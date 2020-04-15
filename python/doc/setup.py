import subprocess

# workaround to install the required LAPACKE library during "Read the Docs" build process
subprocess.run('apt install libblas-dev liblapacke-dev', shell=True, check=True)
