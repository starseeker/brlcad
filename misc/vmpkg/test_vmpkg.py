#!/usr/bin/env python3
#                     T E S T _ V M P K G . P Y
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

"""Offline checks for release inputs, VM commands, and package integrity."""

import argparse
import importlib.util
import io
import json
from pathlib import Path
import stat
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import vmpkg

spec = importlib.util.spec_from_file_location('guest_env', vmpkg.HERE / 'guest-env.py')
guest_env = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guest_env)


class InputsTest(unittest.TestCase):
    def test_kickstart_is_explicit(self):
        rendered = vmpkg.kickstart('ssh-ed25519 YWJj comment')
        self.assertNotIn('@BASE_PACKAGES@', rendered)
        self.assertIn('%packages --nocore --excludedocs --excludeWeakdeps', rendered)
        self.assertIn('sshkey --username=builder "ssh-ed25519 YWJj"', rendered)
        self.assertNotIn('@core', rendered)

    def test_key_cannot_inject_kickstart(self):
        with self.assertRaises(ValueError):
            vmpkg.kickstart('ssh-ed25519 YWJj\n%post')

    def test_boot_menu_preserves_stage2(self):
        rendered = vmpkg.unattended_menu('set default="1"\nset timeout=60\n'
                                         'linuxefi /vmlinuz inst.stage2=hd:LABEL=Rocky quiet\n', True)
        self.assertIn('set default="0"', rendered)
        self.assertIn('inst.stage2=hd:LABEL=Rocky', rendered)
        self.assertIn('inst.ks=hd:LABEL=OEMDRV:/ks.cfg', rendered)

    def test_read_only_boot_menu_is_rewritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'grub.cfg'
            path.write_text('set default="1"\nset timeout=60\n'
                            'linuxefi /vmlinuz inst.stage2=hd:LABEL=Rocky quiet\n')
            read_only_mode = stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH
            path.chmod(read_only_mode)
            vmpkg.rewrite_boot_menu(path, True)
            self.assertIn('inst.ks=hd:LABEL=OEMDRV:/ks.cfg', path.read_text())
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), read_only_mode)

    def test_checked_download_rejects_corruption(self):
        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            (cache / 'input').write_bytes(b'corrupt')
            with self.assertRaises(ValueError):
                vmpkg.download({'url': 'https://example.invalid/input', 'sha256': '0' * 64}, cache)

    def test_vm_commands_use_same_raw_disk_and_seed(self):
        args = argparse.Namespace(work=Path('/tmp/guest'), firmware=Path('/tmp/code.fd'),
                                  cpus=4, memory_mib=24576, install=True,
                                  ssh_port=2222, tap='tap0', name='release-builder')
        qemu = [str(value) for value in vmpkg.qemu_command(args)]
        bhyve = [str(value) for value in vmpkg.bhyve_command(args)]
        self.assertIn('-enable-kvm', qemu)
        self.assertIn('3:0,virtio-blk,/tmp/guest/disk.raw', bhyve)
        self.assertTrue(any('seed.iso' in value for value in qemu))
        self.assertTrue(any('seed.iso' in value for value in bhyve))
        args.install = False
        self.assertFalse(any('seed.iso' in str(value) for value in vmpkg.qemu_command(args)))
        self.assertFalse(any('install.iso' in str(value) for value in vmpkg.bhyve_command(args)))

    def test_archive_escape_rejected_before_extraction(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / 'sources.tar.gz'
            with tarfile.open(archive, 'w:gz') as output:
                member = tarfile.TarInfo('../escape')
                member.size = 1
                output.addfile(member, io.BytesIO(b'x'))
            with self.assertRaises(ValueError):
                vmpkg.extract_archive(archive, root / 'sources')
            self.assertFalse((root / 'escape').exists())

    def test_both_tar_compressions_extract(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for compression in ('gz', 'bz2'):
                archive = root / ('input.tar.' + compression)
                with tarfile.open(archive, 'w:' + compression) as output:
                    member = tarfile.TarInfo('bin/example')
                    member.size = 1
                    output.addfile(member, io.BytesIO(b'x'))
                destination = root / compression
                vmpkg.extract_archive(archive, destination)
                self.assertEqual((destination / 'bin/example').read_bytes(), b'x')

    def test_archive_symlink_escape_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / 'sources.tar.gz'
            with tarfile.open(archive, 'w:gz') as output:
                member = tarfile.TarInfo('link')
                member.type = tarfile.SYMTYPE
                member.linkname = '../../escape'
                output.addfile(member)
            with self.assertRaises(ValueError):
                vmpkg.extract_archive(archive, root / 'sources')


class RpmLockTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.bundle = Path(self.temporary.name)
        key = self.bundle / 'test-key'
        key.write_bytes(b'pinned public key')
        inputs = vmpkg.read_inputs()
        self.inputs = inputs
        inputs['epel_key'] = {'url': 'https://example.invalid/test-key', 'sha256': vmpkg.digest(key)}
        key_policy = patch.object(guest_env, 'read_inputs', return_value=inputs)
        key_policy.start()
        self.addCleanup(key_policy.stop)
        (self.bundle / 'rpms').mkdir()
        rpm = self.bundle / 'rpms/example.rpm'
        rpm.write_bytes(b'pinned RPM bytes')
        self.manifest = {'schema': 1, 'role': 'runtime', 'base_policy': guest_env.base_policy(),
                         'roots': vmpkg.package_roots(vmpkg.HERE / 'runtime-packages.txt'),
                         'repository_config_sha256': vmpkg.digest(vmpkg.HERE / 'rocky.repo'),
                         'rpms': {'example.rpm': {'sha256': vmpkg.digest(rpm), 'nevra': 'example'}}}
        (self.bundle / 'lock.json').write_text(json.dumps(self.manifest))

    def test_valid_bundle(self):
        self.assertEqual(guest_env.validate_bundle(self.bundle), self.manifest)

    def make_builder_bundle(self):
        self.manifest['role'] = 'builder'
        self.manifest['roots'] = vmpkg.package_roots(vmpkg.HERE / 'builder-packages.txt')
        for name in guest_env.BUILDER_ARCHIVES:
            archive = self.bundle / (name + '.archive')
            archive.write_bytes(('pinned ' + name).encode())
            checksum = vmpkg.digest(archive)
            self.inputs[name]['sha256'] = checksum
            self.manifest[name] = {'file': archive.name, 'sha256': checksum}
        (self.bundle / 'lock.json').write_text(json.dumps(self.manifest))

    def test_valid_builder_bundle(self):
        self.make_builder_bundle()
        self.assertEqual(guest_env.validate_bundle(self.bundle), self.manifest)

    def test_modified_ninja_archive_rejected(self):
        self.make_builder_bundle()
        (self.bundle / self.manifest['ninja']['file']).write_bytes(b'modified')
        with self.assertRaisesRegex(ValueError, 'ninja input'):
            guest_env.validate_bundle(self.bundle)

    def test_missing_ninja_pin_rejected(self):
        self.make_builder_bundle()
        del self.manifest['ninja']
        (self.bundle / 'lock.json').write_text(json.dumps(self.manifest))
        with self.assertRaisesRegex(ValueError, 'ninja input'):
            guest_env.validate_bundle(self.bundle)

    def test_builder_archives_cannot_overwrite_each_other(self):
        archives = {}
        for name in guest_env.BUILDER_ARCHIVES:
            path = self.bundle / name / 'archive'
            path.parent.mkdir()
            path.write_bytes(name.encode())
            self.inputs[name]['sha256'] = vmpkg.digest(path)
            archives[name] = path
        args = argparse.Namespace(role='builder', bundle=self.bundle / 'new-bundle', **archives)
        with self.assertRaisesRegex(ValueError, 'basenames must be distinct'):
            guest_env.lock(args)
        self.assertFalse(args.bundle.exists())

    def test_modified_rpm_rejected(self):
        (self.bundle / 'rpms/example.rpm').write_bytes(b'modified')
        with self.assertRaises(ValueError):
            guest_env.validate_bundle(self.bundle)

    def test_modified_key_rejected(self):
        (self.bundle / 'test-key').write_bytes(b'modified key')
        with self.assertRaises(ValueError):
            guest_env.validate_bundle(self.bundle)

    def test_extra_rpm_rejected(self):
        (self.bundle / 'rpms/extra.rpm').write_bytes(b'extra')
        with self.assertRaises(ValueError):
            guest_env.validate_bundle(self.bundle)

    def test_policy_drift_rejected(self):
        with patch.object(guest_env, 'package_roots', return_value=['different']):
            with self.assertRaises(ValueError):
                guest_env.validate_bundle(self.bundle)


if __name__ == '__main__':
    unittest.main()
