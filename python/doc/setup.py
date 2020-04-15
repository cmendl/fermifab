import subprocess

# workaround to install the required LAPACKE library during "Read the Docs" build process
subprocess.run('sudo apt-get install libblas-dev liblapacke-dev', shell=True, check=True)
