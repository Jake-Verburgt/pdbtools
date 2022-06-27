from contextlib import contextmanager
import os
import random
import shutil

@contextmanager
def work_in_dir(dir):
    """Provides a context manager for working in a directory without worrying
       about changing back to the previous wording directory. Meant to be used
       in the context as: with work_in_dir('/path/to/dir'): """
    try:
        old_wd = os.getcwd()
        os.chdir(dir)
        yield
    finally:
        os.chdir(old_wd)


@contextmanager
def in_temp_dir(verbose = False, delete = True):
    """Provides a context manager for working in a directory without worrying
       about deleting it """
    try:
        old_wd = os.getcwd()
        rand = random.randint(10000000, 100000000)
        temp = os.path.join(old_wd, str(rand))
        while os.path.isdir(temp):
            rand = random.randint(10000000, 100000000)
            temp = os.path.join(old_wd, str(rand))
        os.makedirs(temp)
        os.chdir(temp)
        if verbose:
            print(temp)
        yield old_wd
    finally:
        os.chdir(old_wd)
        if delete:
            shutil.rmtree(temp)
            #os.removedirs(temp) #Will not delete directories with stuff in it

