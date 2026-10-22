#!/usr/bin/env python3
#                     V M P K G . P Y
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

"""Prepare release-builder media and run it with bhyve or QEMU/KVM."""

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import select
import shlex
import shutil
import stat
import subprocess
import tarfile
import tempfile
import urllib.request

HERE = Path(__file__).resolve().parent
DEFAULT_CPUS = 4
DEFAULT_MEMORY_MIB = 24576
DEFAULT_DISK_GIB = 200
DEFAULT_SSH_PORT = 2222
BLOCK_BYTES = 1024 * 1024
GIB_BYTES = 1024 ** 3
MIN_DISK_GIB = 16
XVFB_START_SECONDS = 15
XVFB_STOP_SECONDS = 10
XVFB_SCREEN = "800x600x24"


def run(command, dry_run=False, check=True, **kwargs):
    command = [str(arg) for arg in command]
    print('+ ' + ' '.join(shlex.quote(arg) for arg in command), flush=True)
    if not dry_run:
        return subprocess.run(command, check=check, **kwargs)
    return None


def digest(path):
    result = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(BLOCK_BYTES), b''):
            result.update(chunk)
    return result.hexdigest()


def extract_archive(archive, destination):
    with tarfile.open(archive, 'r:*') as source:
        # Git snapshots may contain symlinks. Require link targets to stay
        # within the archive, and reject special files before extraction.
        for member in source.getmembers():
            path = (destination / member.name).resolve()
            if not path.is_relative_to(destination) or member.isdev() or member.isfifo():
                raise ValueError('Unsafe source archive entry: ' + member.name)
            if member.issym() or member.islnk():
                base = path.parent if member.issym() else destination
                if not (base / member.linkname).resolve().is_relative_to(destination):
                    raise ValueError('Source archive link escapes its root: ' + member.name)
        source.extractall(destination, filter='data')


@contextmanager
def virtual_display(log_path):
    previous = os.environ.get('DISPLAY')
    read_fd, write_fd = os.pipe()
    with log_path.open('w') as log:
        try:
            server = subprocess.Popen(['Xvfb', '-displayfd', str(write_fd), '-screen', '0',
                                       XVFB_SCREEN, '-nolisten', 'tcp'],
                                      pass_fds=(write_fd,), stdout=log, stderr=log)
        except OSError:
            os.close(read_fd)
            raise
        finally:
            os.close(write_fd)
        try:
            ready, _, _ = select.select([read_fd], [], [], XVFB_START_SECONDS)
            display = os.read(read_fd, 32).decode().strip() if ready else ''
            if not display.isdigit():
                raise ValueError('Xvfb failed to start; see ' + str(log_path))
            os.environ['DISPLAY'] = ':' + display
            yield
        finally:
            os.close(read_fd)
            server.terminate()
            try:
                server.wait(timeout=XVFB_STOP_SECONDS)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait()
            if previous is None:
                os.environ.pop('DISPLAY', None)
            else:
                os.environ['DISPLAY'] = previous


@contextmanager
def guest_lock(work):
    with (work / '.running.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def external_path(path):
    path = path.resolve()
    if path.is_relative_to(HERE.parents[1]):
        raise ValueError('Store VM media, caches, and source archives outside the source repository')
    return path


def read_inputs():
    return json.loads((HERE / 'inputs.json').read_text())


def download(item, cache):
    target = cache / item['url'].rsplit('/', 1)[1]
    if not target.exists():
        partial = target.with_name(target.name + '.partial')
        with urllib.request.urlopen(item['url'], timeout=120) as source, partial.open('wb') as output:
            shutil.copyfileobj(source, output, BLOCK_BYTES)
        if digest(partial) != item['sha256']:
            partial.unlink()
            raise ValueError('Checksum mismatch for ' + item['url'])
        partial.rename(target)
    if digest(target) != item['sha256']:
        raise ValueError('Cached file checksum mismatch: ' + str(target))
    return target


def fetch(args):
    cache = external_path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    for item in read_inputs().values():
        print(download(item, cache))


def package_roots(path):
    return [line.split('#', 1)[0].strip() for line in path.read_text().splitlines()
            if line.split('#', 1)[0].strip()]


def kickstart(key):
    fields = key.strip().split()
    if len(fields) < 2 or not re.fullmatch(r'(ssh-|ecdsa-|sk-)[A-Za-z0-9@.-]+', fields[0]):
        raise ValueError('An OpenSSH public key is required')
    if not re.fullmatch(r'[A-Za-z0-9+/=]+', fields[1]) or '\n' in key.strip():
        raise ValueError('Invalid public key')
    template = (HERE / 'kickstart.cfg.in').read_text()
    return template.replace('@SSH_KEY@', ' '.join(fields[:2])).replace(
        '@BASE_PACKAGES@', '\n'.join(package_roots(HERE / 'base-packages.txt')))


def unattended_menu(text, grub):
    # Keep the distribution's stage2 label and bootloader layout intact.
    text = re.sub(r'(?m)^(\s*(?:linuxefi|linux|append)\s+.*)$',
                  r'\1 inst.ks=hd:LABEL=OEMDRV:/ks.cfg console=tty0 console=ttyS0,115200n8', text)
    if grub:
        text = re.sub(r'(?m)^set timeout=.*$', 'set timeout=1', text)
        text = re.sub(r'(?m)^set default=.*$', 'set default="0"', text)
    else:
        text = re.sub(r'(?m)^timeout\s+.*$', 'timeout 10', text)
        text = re.sub(r'(?m)^default\s+.*$', 'default linux', text)
    # Avoid the media-test entry, which can add hours to a serial installation.
    text = text.replace('rd.live.check', '')
    return text


def rewrite_boot_menu(path, grub):
    text = unattended_menu(path.read_text(), grub)
    mode = stat.S_IMODE(path.stat().st_mode)
    # ISO extraction preserves read-only permissions. Keep that mode in the
    # remastered image while allowing the host copy to be edited.
    path.chmod(mode | stat.S_IWUSR)
    try:
        path.write_text(text)
    finally:
        path.chmod(mode)


def prepare(args):
    work = external_path(args.work)
    if args.disk_gib < MIN_DISK_GIB:
        raise ValueError('The guest disk must have at least {} GiB virtual capacity'.format(MIN_DISK_GIB))
    if work.exists():
        raise ValueError('Use a new work directory: ' + str(work))
    inputs = read_inputs()
    iso = args.cache.resolve() / inputs['rocky']['url'].rsplit('/', 1)[1]
    if not iso.is_file() or digest(iso) != inputs['rocky']['sha256']:
        raise ValueError('Run fetch first; the pinned Rocky DVD is missing or invalid')
    key = args.public_key.read_text()
    run(['ssh-keygen', '-l', '-f', args.public_key], stdout=subprocess.DEVNULL)
    rendered = kickstart(key)
    work.mkdir(parents=True)
    (work / '.running.lock').touch()
    seed = work / 'seed'
    seed.mkdir()
    shutil.copytree(HERE, seed / 'vmpkg', ignore=shutil.ignore_patterns('__pycache__', '.gitignore'))
    (seed / 'ks.cfg').write_text(rendered)
    run(['xorriso', '-as', 'mkisofs', '-quiet', '-J', '-r', '-V', 'OEMDRV',
         '-o', work / 'seed.iso', seed])
    maps = []
    for entry, grub in [('isolinux/isolinux.cfg', False), ('EFI/BOOT/grub.cfg', True)]:
        target = work / Path(entry).name
        run(['xorriso', '-osirrox', 'on', '-indev', iso, '-extract', '/' + entry, target])
        rewrite_boot_menu(target, grub)
        maps.extend(['-map', target, '/' + entry])
    run(['xorriso', '-indev', iso, '-outdev', work / 'install.iso',
         '-boot_image', 'any', 'replay'] + maps)
    with (work / 'disk.raw').open('xb') as disk:
        disk.truncate(args.disk_gib * GIB_BYTES)
    (work / 'inputs.json').write_text(json.dumps(inputs, indent=2) + '\n')
    print('Prepared installer and sparse disk in ' + str(work))


def firmware_copy(source, destination, dry_run):
    if not destination.exists():
        print('Copy firmware variables: {} -> {}'.format(source, destination))
        if not dry_run:
            shutil.copyfile(source, destination)


def qemu_command(args):
    work = args.work.resolve()
    command = ['qemu-system-x86_64', '-enable-kvm', '-machine', 'q35',
               '-cpu', 'host', '-smp', args.cpus, '-m', args.memory_mib,
               '-nographic', '-no-reboot',
               '-drive', 'if=pflash,format=raw,readonly=on,file=' + str(args.firmware.resolve()),
               '-drive', 'if=pflash,format=raw,file=' + str(work / 'uefi-vars.fd'),
               '-drive', 'file=' + str(work / 'disk.raw') + ',format=raw,if=virtio,discard=unmap',
               '-netdev', 'user,id=net0,hostfwd=tcp:127.0.0.1:{}-:22'.format(args.ssh_port),
               '-device', 'virtio-net-pci,netdev=net0']
    if args.install:
        command += ['-drive', 'file=' + str(work / 'install.iso') + ',media=cdrom,readonly=on',
                    '-drive', 'file=' + str(work / 'seed.iso') + ',media=cdrom,readonly=on',
                    '-boot', 'once=d']
    return command


def bhyve_command(args):
    work = args.work.resolve()
    command = ['bhyve', '-AHPD', '-c', args.cpus, '-m', str(args.memory_mib) + 'M',
               '-s', '0:0,hostbridge', '-s', '1:0,lpc', '-s', '2:0,virtio-net,' + args.tap,
               '-s', '3:0,virtio-blk,' + str(work / 'disk.raw'),
               '-l', 'com1,stdio', '-l', 'bootrom,' + str(args.firmware.resolve()) + ',' + str(work / 'uefi-vars.fd')]
    if args.install:
        command += ['-s', '4:0,ahci-cd,' + str(work / 'install.iso'),
                    '-s', '5:0,ahci-cd,' + str(work / 'seed.iso')]
    return command + [args.name]


def launch(args):
    expected_host = 'Linux' if args.action == 'qemu' else 'FreeBSD'
    if not args.dry_run and platform.system() != expected_host:
        raise ValueError(args.action + ' requires a ' + expected_host + ' host')
    for value in (args.cpus, args.memory_mib):
        if value <= 0:
            raise ValueError('CPU count and memory must be positive')
    # Commas delimit QEMU and bhyve device options, even without a shell.
    for path in (args.work, args.firmware, args.vars_template):
        if ',' in str(path):
            raise ValueError('VM device paths cannot contain commas')
    if not args.dry_run and not (args.work / 'disk.raw').is_file():
        raise ValueError('Run prepare first')
    command = qemu_command(args) if args.action == 'qemu' else bhyve_command(args)
    if args.dry_run:
        run(command, dry_run=True)
        return
    vm_device = None
    if args.action == 'bhyve':
        if not re.fullmatch(r'[A-Za-z0-9_-]+', args.name) or ',' in args.tap:
            raise ValueError('Invalid bhyve guest name or tap interface')
        vm_device = Path('/dev/vmm') / args.name
        if vm_device.exists():
            raise ValueError('A bhyve guest with this name already exists')
    elif not 0 < args.ssh_port < 65536:
        raise ValueError('Invalid SSH forwarding port')
    with guest_lock(args.work):
        firmware_copy(args.vars_template, args.work / 'uefi-vars.fd', False)
        try:
            result = run(command, check=False)
            allowed_status = (0, 1, 2) if args.action == 'bhyve' else (0,)
            if result.returncode not in allowed_status:
                raise subprocess.CalledProcessError(result.returncode, command)
        finally:
            if vm_device is not None and vm_device.exists():
                run(['bhyvectl', '--destroy', '--vm=' + args.name])


def clone(args):
    destination = external_path(args.output)
    if destination.exists():
        raise ValueError('Use a new clone directory')
    with guest_lock(args.work):
        destination.mkdir(parents=True)
        (destination / '.running.lock').touch()
        run(['qemu-img', 'convert', '-f', 'raw', '-O', 'raw', '-S', '4k',
             args.work / 'disk.raw', destination / 'disk.raw'])
        for name in ('uefi-vars.fd', 'inputs.json'):
            shutil.copyfile(args.work / name, destination / name)


def snapshot_one(source, target, origin):
    source = source.resolve()
    run(['git', 'clone', '--no-local', '--depth', '1', source, target], stdout=subprocess.DEVNULL)
    run(['git', '-C', target, 'remote', 'set-url', 'origin', origin])
    result = subprocess.check_output(['git', '-C', source, 'ls-files', '--cached', '--others',
                                      '--exclude-standard', '-z'])
    for raw in result.split(b'\0'):
        if not raw:
            continue
        name = os.fsdecode(raw)
        original, copied = source / name, target / name
        if original.is_dir():
            # Gitlinks are populated by bext on demand in the Linux guest.
            continue
        if copied.exists() or copied.is_symlink():
            copied.unlink()
        if original.exists() or original.is_symlink():
            copied.parent.mkdir(parents=True, exist_ok=True)
            if original.is_symlink():
                copied.symlink_to(os.readlink(original))
            else:
                shutil.copy2(original, copied)
    return subprocess.check_output(['git', '-C', target, 'rev-parse', 'HEAD']).decode().strip()


def snapshot(args):
    args.output = external_path(args.output)
    if args.output.exists():
        raise ValueError('Source archive already exists')
    with tempfile.TemporaryDirectory(prefix='vmpkg-sources-') as temporary:
        root = Path(temporary)
        revisions = {}
        for name, source in [('brlcad', args.brlcad), ('bext', args.bext)]:
            revisions[name] = snapshot_one(source, root / name, 'https://github.com/BRL-CAD/' + name)
        (root / 'revisions.json').write_text(json.dumps(revisions, indent=2) + '\n')
        with tarfile.open(args.output, 'w:gz') as archive:
            for name in ('brlcad', 'bext', 'revisions.json'):
                archive.add(root / name, arcname=name)
    print(digest(args.output) + '  ' + str(args.output))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest='action', required=True)
    fetch_parser = actions.add_parser('fetch', help='Retrieve and verify pinned ISO and CMake inputs')
    fetch_parser.add_argument('--cache', type=Path, required=True)
    fetch_parser.set_defaults(function=fetch)
    prepare_parser = actions.add_parser('prepare', help='Create unattended installer media and a fresh disk')
    prepare_parser.add_argument('--cache', type=Path, required=True)
    prepare_parser.add_argument('--work', type=Path, required=True)
    prepare_parser.add_argument('--public-key', type=Path, required=True)
    prepare_parser.add_argument('--disk-gib', type=int, default=DEFAULT_DISK_GIB)
    prepare_parser.set_defaults(function=prepare)
    for name in ('qemu', 'bhyve'):
        launch_parser = actions.add_parser(name, help='Run the installer or the installed guest')
        launch_parser.add_argument('--work', type=Path, required=True)
        launch_parser.add_argument('--firmware', type=Path, required=True)
        launch_parser.add_argument('--vars-template', type=Path, required=True)
        launch_parser.add_argument('--cpus', type=int, default=DEFAULT_CPUS)
        launch_parser.add_argument('--memory-mib', type=int, default=DEFAULT_MEMORY_MIB)
        launch_parser.add_argument('--install', action='store_true')
        launch_parser.add_argument('--dry-run', action='store_true')
        if name == 'qemu':
            launch_parser.add_argument('--ssh-port', type=int, default=DEFAULT_SSH_PORT)
        else:
            launch_parser.add_argument('--tap', required=True)
            launch_parser.add_argument('--name', default='release-builder')
        launch_parser.set_defaults(function=launch)
    clone_parser = actions.add_parser('clone', help='Sparsely clone a shut-down base or prepared guest')
    clone_parser.add_argument('--work', type=Path, required=True)
    clone_parser.add_argument('--output', type=Path, required=True)
    clone_parser.set_defaults(function=clone)
    snapshot_parser = actions.add_parser('snapshot', help='Export current sources, including local changes')
    snapshot_parser.add_argument('--brlcad', type=Path, required=True)
    snapshot_parser.add_argument('--bext', type=Path, required=True)
    snapshot_parser.add_argument('--output', type=Path, required=True)
    snapshot_parser.set_defaults(function=snapshot)
    args = parser.parse_args()
    try:
        args.function(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, str(error) + '\n')


if __name__ == '__main__':
    main()
