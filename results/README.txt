Instructions to Run the quantum-open-system-polynomial-sid reproducible Docker Image

1. Download the Docker Image:
1.1 Create a separate directory to keep things organized:

    mkdir oqsid
    cd oqsid

1.2 Download the oqsid Docker image from the following link:

https://drive.google.com/file/d/1a61hLOgV6hq4TQK3LhxJmHGQ8w0s0FYq/view?usp=sharing

1.3 Load the Docker Image from the .tar File

    docker load -i oqsid.tar

2. Obtain a Free MOSEK License:
2.1 Visit MOSEK's Academic License Page to obtain a license:
https://www.mosek.com/products/academic-licenses/

2.2 Save the mosek.lic file to the same directory where the Docker image is located.

3. Run the Docker Image:
3.1 Navigate to the directory containing the mosek.lic file.
 Then, mount this directory as the /root/OQSID/results directory in the container and run the Docker image:

    docker run -v "$(pwd)":/root/OQSID/results oqsid

4. What Happens Next:
The Docker image will automatically:

    Download the source code and input data from the repository:
    https://github.com/zpopovych/OQSID_for_Docker
    Run the necessary computations.
    Save the resulting data files and figures to your current directory.

Note:
Ensure that Docker is installed and running on your system before proceeding. 
For more information, visit the Docker installation page.
https://docs.docker.com/get-started/get-docker/