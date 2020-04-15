import subprocess

# workaround to install the required LAPACKE library during "Read the Docs" build process
subprocess.run('apt-get -y install sudo', shell=True, check=True)
subprocess.run('sudo apt-get -y install libblas-dev liblapacke-dev', shell=True, check=True)
