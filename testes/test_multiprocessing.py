from multiprocessing import Process, current_process, Array, Pool, Queue, Manager
import threading
from concurrent.futures import ProcessPoolExecutor

import multiprocessing

def f(name):
    print("hello", name)


if __name__ == "__main__":
    p = Process(target=f, args=("bob",))

#print("Number of cpu : ", multiprocessing.cpu_count())



# Pool():  Determina os processos que vão pra cada núcleo de computacão disponíveis
#    starmap()  [(1,2), (3, 4)] results in [func(1,2), func(3,4)]
#    apply(func, arg=())  uses just one process
# Process():  Envia os processos sequencialmente para os núcleos
# start()  inicia o processo, roda, cria e inicia a thread principal
#          a thread principal pode ficar bloqueada por outros processos
#
# daemon: background process This means that daemon child processes can run in the
#         background and do not prevent a Python program from exiting when the “main
#         parts of the program” have finished
# process.join() vai esperar o processo terminar

# Process: An instance of the Python interpreter has at least one thread called the MainThread.
# Thread: A thread of execution within a Python process, such as the MainThread or a new thread.

# Mutex: mutual exclusion (mutex) lock in Python via the multiprocessing.Lock class.
#        mutex lock is used to protect critical sections of code from concurrent execution
#        synchronization primitive intended to prevent a race condition
#        from multiprocessing import Lock

# Array: A shared C-type memory to be used for the processes
from multiprocessing import Process
import math

def task3(arg):
    x = sum([math.sqrt(i) for i in range(1, arg)])
    return x

# protect the entry point
if __name__ == '__main__':
    # report a message
    print('Starting task...')
    # create the process pool
    with ProcessPoolExecutor(4) as exe:
        # perform calculations
        results = exe.map(task3, range(1,5000))
        # Como pega o resultado?????
    # report a message
    print('Done.')

# define a cpu-intensive task
def task2(arg):
    return sum([math.sqrt(i) for i in range(1, arg)])
# protect the entry point
if __name__ == '__main__':
    # report a message
    print('Starting task...')
    # create the process pool with 4 
    with Pool(4) as pool:
        # perform calculations
        results = pool.map(task2, range(1,5000))
        #pool.close()?
    # report a message
    print('Done.')  # This is faster than ProcessPoolExecutor


# Pesquisar: multiprocessing pool methods
"""
# SuperFastPython.com
# list all active child processes
from time import sleep
from multiprocessing import active_children
# function to execute in a new process
def task1(i):
    # block for a moment
    if i ==2:
        sleep(15)
    else:
        sleep(3)
    print(f'Process {i} done', flush=True)

# entry point
if __name__ == '__main__':
    # create a number of child processes
    processes = [Process(target=task1, args=(i,)) for i in range(5)]
    # start the child processes
    for process in processes:
        process.start()
    # get a list of all active child processes
    children = active_children()
    for process in processes:
        process.join()
    # report a count of active children
    print(f'Active Children Count: {len(children)}')
    # report each in turn
    for child in children:
        print(child)


"""

x = list(range(5,10))
def f2(d, l, i, x):
    d[i] = i
    x[i] = i
    #d['2'] = 2
    #d[0.25] = None
    l.reverse()
    if i == 2:
        sleep(15)
    else:
        sleep(5)
    print(d, i, flush=True)

if __name__ == '__main__':
    with Manager() as manager:
        d = manager.dict()
        l = manager.list(range(10))
        
        for i in range(5):
            p = Process(target=f2, args=(d, l, i, x))
            p.start()
        p.join()
        print(d)
        print(l)


        
# entry point
if __name__ == "__main__":
    
    # enable support for multiprocessing
    multiprocessing.freeze_support()
    # create the process
    process = Process()
    # report the process identifier
    print(process.pid)
    print(process.is_alive())
    # start the process
    process.start()
    print(process.is_alive())
    # report the process identifier
    print(process.pid)
    print(process.exitcode)
    process.join()
    print("2", process.is_alive())
    print(process.exitcode)
    print(process.name)
    
    


# example of using a barrier with processes
from time import sleep
from random import random
from multiprocessing import Process
from multiprocessing import Barrier
# target function to prepare some work
def task(barrier, number, shared):
    
    # generate a unique value
    value = random()
    # block for a moment
    sleep(value)
    shared[number] *= 2
    # report result
    print(f'Process {number} done, got: {value}', flush=True)
    # wait on all other processes to complete
    barrier.wait()
# entry point
if __name__ == '__main__':
    
    # enable support for multiprocessing
    multiprocessing.freeze_support()   
    # create a barrier
    barrier = Barrier(5 + 1)
    
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    shared = Array('f', a)
    # create the worker processes
    for i in range(5):
        # start a new process to perform some work
        worker = Process(target=task, args=(barrier, i, shared))
        worker.start()
    # wait for all processes to finish
    print('Main process waiting on all results...')
    barrier.wait()

    print(shared)
    # report once all processes are done
    print('All processes have their result')


def f1(q):
    q.put([42, None, 'hello'])

if __name__ == '__main__':

    q = Queue()
    for i in range(5):
        p = Process(target=f1, args=(q,))
        p.start()
    for i in range(5):
        p.join()
    
    print(q.get()) # prints "[42, None, 'hello']"
    print(q.get()) # prints "[42, None, 'hello']"5
    #with Pool(15) as p:
    #    print(p.map(f, [1, 2, 3, 4, 5, 6]))


