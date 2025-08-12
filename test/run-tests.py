#!/usr/bin/env python3

import argparse
import configparser
import os
import sys
from pathlib import Path
import shutil
import glob
import subprocess
from subprocess import CalledProcessError
import importlib.util
from collections import OrderedDict
import datetime

pluto_dir = os.environ['PLUTO_DIR']
sys.path.append(pluto_dir + '/Tools/Python/')
from make_problem import MakeProblem

debug = False

class Test:
    def __init__(self, directory, makefile_defs):
        if os.path.isdir(directory):
            self.directory = Path(directory)
            self.config_file = self.directory / 'test.ini'
        elif os.path.isfile(directory):
            self.config_file = directory
            self.directory = Path(os.path.dirname(directory))

        self.name = Path(self.directory.stem) / self.config_file.stem
        self.tmp_dir = self.directory / 'test'

        config = configparser.ConfigParser()
        config.read(self.config_file)

        self.definitions = config['Definitions']
        self.definitions_user = config['DefinitionsUser']

        self.physics = self.definitions['PHYSICS']
        self.geometry = self.definitions['GEOMETRY']
        self.dimensions = self.definitions['DIMENSIONS']
        self.body_force = self.definitions['BODY_FORCE']

        self.long_name = config['Info']['Name']
        self.description = config['Info']['Description']

        self.extra_files = []
        if 'ExtraFiles' in config['Info']:
            self.extra_files = [self.directory / x for x in config['Info']['ExtraFiles'].split(',')]

        self.pluto_arguments = []
        if 'PlutoArguments' in config['Info']:
            self.pluto_arguments = config['Info']['PlutoArguments'].strip().split(' ')

        self.verify_script = self.directory / config['Verify']['Script']

        self.makefile_defs = makefile_defs

    def __str__(self):
        return f'{self.name}'

    def print_test_info(self):
        print(f'[{self.physics}, {self.dimensions}D, {self.geometry}, {self.body_force}] {self.name}: {self.long_name}\n{self.description}')

def setup_test(cur_test):
    print(f'Setting up {cur_test}')

    # clean up any old test data if it still exists
    if cur_test.tmp_dir.exists():
        cleanup_test(cur_test)

    cur_test.tmp_dir.mkdir() # create test directory

    # copy source files
    shutil.copytree("../src/", cur_test.tmp_dir, dirs_exist_ok=True)
    # for f in glob.glob("../src/*"):
    #     print(f'Copying {f}')
    #     shutil.copy(f, cur_test.tmp_dir)

    shutil.copy(cur_test.directory / 'pluto.ini', cur_test.tmp_dir)
    shutil.copy(cur_test.directory / 'makefile', cur_test.tmp_dir)

    # copy extra files
    for f in cur_test.extra_files:
        shutil.copy(f, cur_test.tmp_dir)

    # edit source files

    # definitions.h
    with open(cur_test.tmp_dir / 'definitions.h', 'r') as file:
        filedata = file.readlines()

    for v in cur_test.definitions:
        line_number = next(i for i, string in enumerate(filedata) if v.upper() in string)
        filedata[line_number] = '#define  ' + v.upper().ljust(28) + '   ' + f'{cur_test.definitions[v].upper()}\n'

    with open(cur_test.tmp_dir / 'definitions.h', 'w') as file:
        file.writelines(filedata)

    # definitions_usr.h
    with open(cur_test.tmp_dir / 'definitions_usr.h', 'r') as file:
        filedata = file.readlines()

    for v in cur_test.definitions_user:
        line_number = next(i for i, string in enumerate(filedata) if v.upper() in string)
        filedata[line_number] = '#define  ' + v.upper().ljust(21) + '   ' + f'{cur_test.definitions_user[v].upper()}\n'

    with open(cur_test.tmp_dir / 'definitions_usr.h', 'w') as file:
        file.writelines(filedata)

    # makefile
    if cur_test.makefile_defs is not None:
        with open(cur_test.tmp_dir / 'makefile', 'r') as file:
            filedata = file.readlines()

        line_number = next(i for i, string in enumerate(filedata) if 'ARCH' in string)
        filedata[line_number] = 'ARCH'.ljust(12) + ' = ' + f'{cur_test.makefile_defs}\n'

        with open(cur_test.tmp_dir / 'makefile', 'w') as file:
            file.writelines(filedata)

    # generate our makefile
    MakeProblem(str(cur_test.tmp_dir), pluto_dir, 1, 1)

    # build pluto
    proc = subprocess.run(['make', '-j'], check = True, cwd = cur_test.tmp_dir,  stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

    if debug:
        print(proc.stdout.decode('UTF-8'))

def run_test(cur_test):
    print(f'Running {cur_test}')

    # run pluto
    start_time = datetime.datetime.now()
    proc = subprocess.run(['./pluto'] + cur_test.pluto_arguments, check = True, cwd = cur_test.tmp_dir, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

    end_time = datetime.datetime.now()

    print(f'Test run in {end_time - start_time}')

    if debug:
        print(proc.stdout.decode('UTF-8'))

def verify_test(cur_test):
    # verify test
    spec = importlib.util.spec_from_file_location("test", cur_test.verify_script)
    verify = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verify)

    print('Results:')
    verify.verify(cur_test.tmp_dir)

def cleanup_test(cur_test):
    print(f'Cleaning up {cur_test}')
    if cur_test.tmp_dir.exists():
        shutil.rmtree(cur_test.tmp_dir, ignore_errors = False)

def execute_test(cur_test, *, cleanup = True):
    try:
        cur_test.print_test_info()
        setup_test(cur_test)
        run_test(cur_test)
        verify_test(cur_test)
    except CalledProcessError as e:
        print(f'Command {e.cmd} returned non-zero exit status {e.returncode}')
        if debug:
            print(e.stdout.decode('UTF-8'))
    except Exception as e:
        print(e)
    finally:
        if cleanup:
            cleanup_test(cur_test)
        print('')

def main():
    global debug

    ap = argparse.ArgumentParser()
    ap.add_argument('test_name', nargs = '*')
    ap.add_argument('-n', '--no_cleanup', help = "Don't clean up after a test", action='store_true')
    ap.add_argument('-c', '--clean', help = 'Clean all tests', action='store_true')
    ap.add_argument('-v', '--verify', help = 'Verify tests', action='store_true')
    ap.add_argument('-d', '--debug', help = 'Enable debug mode (more output)', action='store_true')
    ap.add_argument('-l', '--list', help = 'List all available tests', action = 'store_true')
    ap.add_argument('-m', '--makefile_defs', help = 'Select a specific makefile definition')
    args = ap.parse_args()

    debug = args.debug

    base_directory = Path ('.')

    test_files = {}

    for root,dirs,files in os.walk(base_directory):
        for f in files:
            if f.startswith('test') and f.endswith('ini'):
                test_files[str((Path(root) / os.path.splitext(f)[0]).relative_to('.'))] = Path(os.path.join(root, f))

    # sort test entries
    test_files = OrderedDict(sorted(test_files.items()))

    if args.test_name == []:
        tests_to_execute = test_files.keys()
    else:
        tests_to_execute = args.test_name

    makefile_defs = args.makefile_defs

    for t in tests_to_execute:
        tdirs = []
        # if not 'test' in t:
        for k in test_files.keys():
            if t in k:
                tdirs.append(test_files[k])

        for tdir in tdirs:
            cur_test = Test(
                directory = tdir,
                makefile_defs = makefile_defs
            )

            if args.clean:
                cleanup_test(cur_test)
                continue
            if args.verify:
                verify_test(cur_test)
                continue
            if args.list:
                cur_test.print_test_info()
                print("")
                continue

            execute_test(cur_test, cleanup = not args.no_cleanup)

if __name__ == "__main__":
    main()
