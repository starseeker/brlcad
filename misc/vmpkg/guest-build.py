#!/usr/bin/env python3
#                     G U E S T - B U I L D . P Y
# BRL-CAD
#
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided
# with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
# products derived from this software without specific prior written
# permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
###

"""Build a source snapshot and retain packages, provenance, and failure logs."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess

from vmpkg import HERE, digest, extract_archive, virtual_display


def logged(command, name, logs, cwd=None):
    print('+ ' + ' '.join(str(value) for value in command), flush=True)
    with (logs / (name + '.log')).open('a') as output:
        subprocess.run(['/usr/bin/time', '-v', '-o', str(logs / (name + '.resources'))] +
                       [str(value) for value in command], cwd=cwd, check=True,
                       stdout=output, stderr=subprocess.STDOUT)



def build(args):
    if os.geteuid() == 0:
        raise ValueError('Run builds as the builder account')
    spec = importlib.util.spec_from_file_location('guest_env', HERE / 'guest-env.py')
    environment = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(environment)
    receipt = environment.verify()
    if receipt['role'] != 'builder':
        raise ValueError('A builder RPM bundle must be installed before building')
    root = args.work.resolve()
    root.mkdir(parents=True, exist_ok=True)
    fingerprint = digest(args.sources)
    recipe = {'sources_sha256': fingerprint, 'jobs': args.jobs, 'profile': args.profile,
              'lock_sha256': receipt['lock_sha256'],
              'cflags': os.environ.get('CFLAGS', ''), 'cxxflags': os.environ.get('CXXFLAGS', '')}
    identity = root / 'recipe.json'
    if identity.exists():
        if json.loads(identity.read_text()) != recipe:
            raise ValueError('Build inputs or options changed; use a new work directory')
    else:
        if any(root.iterdir()):
            raise ValueError('Use an empty work directory')
        extract_archive(args.sources, root / 'sources')
        identity.write_text(json.dumps(recipe, sort_keys=True, indent=2) + '\n')
    sources = root / 'sources'
    logs = root / 'logs'
    logs.mkdir(exist_ok=True)
    output = root / 'output'
    output.mkdir(exist_ok=True)
    shutil.copyfile(HERE / 'environment-lock.json', output / 'environment-lock.json')
    (output / 'provenance.json').write_text(json.dumps({
        'sources_sha256': fingerprint, 'environment': receipt,
        'build_jobs': args.jobs, 'profile': args.profile,
        'brlcad_ref': subprocess.check_output(['git', '-C', sources / 'brlcad', 'rev-parse', 'HEAD']).decode().strip(),
        'bext_ref': subprocess.check_output(['git', '-C', sources / 'bext', 'rev-parse', 'HEAD']).decode().strip()
    }, indent=2) + '\n')
    if not (sources / 'bext/CMake/BextReleaseBuild.cmake').is_file():
        raise ValueError('Source snapshot needs bext release-build support')
    bext_build, bext_output = root / 'bext-build', root / 'bext-output'
    configure = ['cmake', '-S', sources / 'bext', '-B', bext_build, '-G', 'Ninja',
                 '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_INSTALL_PREFIX=' + str(bext_output),
                 '-DENABLE_ALL=ON', '-DGIT_SHALLOW_CLONE=ON',
                 '-DEXT_BUILD_ALWAYS=OFF', '-DBEXT_DISCARD_WORK=ON',
                 '-DBEXT_BUILD_JOBS=' + str(args.jobs), '-DBEXT_LINK_JOBS=1',
                 '-DPython3_EXECUTABLE=/usr/bin/python3.11']
    if args.profile == 'everything':
        configure += ['-DUSE_APPLESEED=ON', '-DUSE_OSPRAY=ON']
    logged(configure, 'bext-configure', logs)
    logged(['cmake', '--build', bext_build, '--parallel', '1'], 'bext-build', logs)
    brlcad_build = root / 'brlcad-build'
    logged(['cmake', '-S', sources / 'brlcad', '-B', brlcad_build, '-G', 'Ninja',
            '-DCMAKE_BUILD_TYPE=Release', '-DBRLCAD_EXT_DIR=' + str(bext_output),
            '-DBRLCAD_ENABLE_QT=ON', '-DBRLCAD_X11_PROVIDER=SYSTEM',
            '-DPython3_EXECUTABLE=/usr/bin/python3.11'], 'brlcad-configure', logs)
    logged(['cmake', '--build', brlcad_build, '--parallel', str(args.jobs)], 'brlcad-build', logs)
    logged(['cmake', '--build', brlcad_build, '--target', 'regress-repository',
            'regress-licenses', '--parallel', str(args.jobs)], 'repository-checks', logs)
    os.environ['LIBGL_ALWAYS_SOFTWARE'] = '1'
    with virtual_display(logs / 'Xvfb.log'):
        logged(['ctest', '--test-dir', brlcad_build, '--output-on-failure',
                '--parallel', str(args.jobs)], 'tests', logs)
    for generator, suffix in [('RPM', '.rpm'), ('TGZ', '.tar.gz'), ('TBZ2', '.tar.bz2')]:
        staging = root / ('package-' + generator.lower())
        logged(['cpack', '--config', brlcad_build / 'CPackConfig.cmake',
                '-G', generator, '-B', staging], 'package-' + generator.lower(), logs)
        packages = list(staging.glob('*' + suffix))
        if not packages:
            raise ValueError('CPack produced no ' + generator + ' package')
        for package in packages:
            package.replace(output / package.name)
        shutil.rmtree(staging)
    with (output / 'SHA256SUMS').open('w') as checksums:
        for artifact in sorted(output.iterdir()):
            if artifact.name != 'SHA256SUMS' and artifact.is_file():
                checksums.write(digest(artifact) + '  ' + artifact.name + '\n')
    print('Packages: ' + str(output))
    print('Logs: ' + str(logs))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--profile', choices=['bundled', 'everything'], default='bundled')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    try:
        build(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, str(error) + '\nSee the work directory logs for command output.\n')


if __name__ == '__main__':
    main()
