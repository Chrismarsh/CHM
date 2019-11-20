#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, shutil, platform, re
from conans import ConanFile, AutoToolsBuildEnvironment, tools

class GmpConan(ConanFile):
    """ Building GMP for the intention of using it to build CGAL """

    name        = 'gmp'
    version     = '6.1.2'
    md5_hash    = '8ddbb26dc3bd4e2302984debba1406a5'
    description = 'The GNU Multiple Precision Arithmetic Library'
    url         = 'https://github.com/CGAL/conan-gmp'
    license     = 'MIT' # TODO: fix this
    settings    = 'os', 'compiler', 'arch', 'build_type'
    options = {
        'shared':            [True, False],
        'static':            [True, False],
        'disable_assembly':  [True, False],
        'enable_fat':        [True, False],
        'enable_cxx':        [True, False],
        'disable-fft':       [True, False],
        'enable-assert':     [True, False],
        'msvc':              [12, 15],
    }
    default_options = (
        'shared=True',
        'static=True',
        'disable_assembly=False',
        'enable_fat=False',
        'enable_cxx=True',
        'disable-fft=False',
        'enable-assert=False',
        'msvc=12'
    )

    def source(self):
        zip_name = 'gmp-{version}.tar.bz2'.format(version=self.version)
        tools.download('http://gnu.uberglobalmirror.com/gmp/{zip_name}'.format(zip_name=zip_name), zip_name)
        # Alternative
        # tools.download(f'http://gmplib.org/download/gmp/{zip_name}', zip_name)
        tools.check_md5(zip_name, self.md5_hash)
        tools.unzip(zip_name)
        shutil.move('gmp-{version}'.format(version=self.version), 'gmp')
        os.unlink(zip_name)

    def configure(self):
        if 'Windows' == self.settings.os:
            # On Windows, we can only build the static OR shared
            self.options.static = not self.options.shared

    def build(self):
        with tools.chdir(self.name):
            autotools = AutoToolsBuildEnvironment(self, win_bash=(platform.system() == "Windows"))

            env_vars = {}
            args = []

            if 'gcc' == self.settings.compiler and 'Windows' == platform.system():
                args.append('--prefix=%s'%tools.unix_path(self.package_folder))
            else:
                args.append('--prefix=%s'%self.package_folder)

            for option_name in self.options.values.fields:
                activated = getattr(self.options, option_name)
                if not re.match(r'enable|disable', option_name):
                    continue
                if activated:
                    option_name = option_name.replace("_", "-")
                    self.output.info("Activated option! %s"%option_name)
                    args.append('--%s'%option_name)

            args.append('--%s-shared'%('enable' if self.options.shared else 'disable'))
            args.append('--%s-static'%('enable' if self.options.static else 'disable'))

            if self.settings.os == "Linux" or self.settings.os == "Macos":
                autotools.fpic = True
                if self.settings.arch == 'x86':
                    env_vars['ABI'] = '32'
                    autotools.cxx_flags.append('-m32')

            # Debug
            self.output.info('Configure arguments: %s'%' '.join(args))

            # Set up our build environment
            with tools.environment_append(env_vars):
                autotools.configure(args=args)

            autotools.make()
            if 'gcc' == self.settings.compiler and 'Windows' == self.settings.os:
                if self.env['DLLTOOL'] is not None:
                    self.run('{dlltool} --output-lib gmp.lib --input-def .libs/libgmp-3.dll.def --dllname libgmp-10.dll'.format(dlltool=self.env['DLLTOOL']))
            autotools.make(args=['install'])

    def package(self):
        self.copy("COPYING*", src="gmp", dst="")
        self.copy("gmp.lib",  src="gmp", dst="lib")

    def package_info(self):
        # We get lib<lib>.a: mpn, mpz, mpq, mpf, printf, scanf, random,
        #  cxx, gmp, gmpxx
        # Not sure if these should all be added?  I think CGAL only wants
        # the first one

        self.cpp_info.libs = ['gmp']

    def package_id(self):
        # On windows, we cross compile this with mingw.. But because it's
        # compatible with MSVC, set it's hash to reflect that.
        if 'gcc' == self.settings.compiler and 'Windows' == self.settings.os:
            self.info.settings.compiler = 'Visual Studio'
            self.info.settings.compiler.version = int(str(self.options.msvc))

            runtime = 'MD' if self.options.shared else 'MT'
            if self.settings.build_type == 'Debug':
                runtime += 'd'
            self.info.settings.compiler.runtime = runtime

# vim: ts=4 sw=4 expandtab ffs=unix ft=python foldmethod=marker :
