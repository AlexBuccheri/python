from parse_args import run_executable

class Args():
    def __init__(self, omp_num_threads, exe, np, build_type):
        self.omp_num_threads = omp_num_threads
        self.exe = exe
        self.np = np
        self.build_type = [build_type]



args = Args(2, 'test', 2, 'mpi')

output = run_executable(args, '')

print("stdout:")
for line in output['stdout']:
    print(line)

print("stderr:")
for line in output['stderr']:
    print(line)
