# Instructions to Run the Quantum-Open-System-Polynomial-SID Docker Image

## 1. Download the Docker Image

### 1.1 Create a Directory for the Project
Run the following commands to organize your files:
```bash
mkdir oqsid
cd oqsid
```

### 1.2 Download the Docker Image
Access the image via the link below and download it to your created directory:  
[Download OQSID Docker Image](https://drive.google.com/file/d/1a61hLOgV6hq4TQK3LhxJmHGQ8w0s0FYq/view?usp=sharing)

### 1.3 Load the Docker Image
Once downloaded, load the image into Docker using:
```bash
docker load -i oqsid.tar
```

---

## 2. Obtain a Free MOSEK License

### 2.1 Apply for a MOSEK Academic License
Visit MOSEK's academic license page and follow the instructions to get your free license:  
[MOSEK Academic Licenses](https://www.mosek.com/products/academic-licenses/)

### 2.2 Save the License File
Save the `mosek.lic` file to the directory where the Docker image (`oqsid.tar`) is located (e.g., `oqsid`).

---

## 3. Run the Docker Image

### 3.1 Navigate to the Working Directory
Ensure you are in the directory containing the `mosek.lic` file.

### 3.2 Run the Docker Image
Use the following command to mount your current directory as the results folder and start the container:
```bash
docker run -v "$(pwd)":/root/OQSID/results oqsid
```

---

## 4. What Happens Next?
The Docker image will automatically:
- **Download the source code and input data** from the GitHub repository:  
  [OQSID Source Code](https://github.com/zpopovych/OQSID_for_Docker)
- **Perform all necessary computations.**
- **Save results (data files and figures)** to your current working directory.

---

## Additional Notes

- **Ensure Docker is Installed and Running**  
  If Docker isn’t installed on your system, follow the installation guide:  
  [Docker Installation Guide](https://docs.docker.com/get-started/get-docker/)

- **Supported Environments**  
  The Docker image is designed to run seamlessly on Linux, macOS, and Windows systems with Docker installed.

For assistance or troubleshooting, refer to the repository's documentation or reach out via the GitHub Issues page.
