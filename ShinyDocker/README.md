### **Running the Docker Container Using Docker Desktop**:

To embark on Building Docker for BATCH-FLEX-Shiny, we encourage you to install Docker Desktop on your system. Once you've installed Docker Desktop, you can build the Docker image and run the container with the following simple steps:

#### **Shiny Visualization Applications Docker Setup:**

**1. Install Docker Desktop**: Download and install Docker Desktop for your operating system from the official Docker website ([https://www.docker.com/products/docker-desktop).](https://www.docker.com/products/docker-desktop).)

**2. Build the Docker Image**:

-   Clone the Docker project repository from its source.

    Click github.dev

    ![](https://github.com/shawlab-moffitt/BATCH-FLEX-ShinyApp/blob/main/ShinyDocker/BATCH-FLEX-Shiny%20Docker%20Desktop%20images_04.png?raw=true)

    Download ShinyDocker folder

![](https://github.com/shawlab-moffitt/BATCH-FLEX-ShinyApp/blob/main/ShinyDocker/BATCH-FLEX-Shiny%20Docker%20Desktop%20images_05.png?raw=true)

-   Open Windows command prompt and navigate to the directory containing the Dockerfie.

```         
C:\>cd C:\Users\Administrator\Desktop\ShinyDocker
```

-   Use the `docker build` command to build the Docker image. For example:

```         
docker build -t batch_flex_shinyapp .
```

**3. Run the Docker Container**:

-   Once the Docker image is successfully built, you can run the image using Docker Desktop. For instance, click 3838:3838 to display web page as following.

    ![](https://github.com/chingyaousf/Introduction-to-Building-Docker-for-BERLIN-Pipeline/blob/main/data/BATCH-FLEX-Shiny%20Docker%20Desktop%20images_01.png?raw=true)

    ![](https://github.com/chingyaousf/Introduction-to-Building-Docker-for-BERLIN-Pipeline/blob/main/data/BATCH-FLEX-Shiny%20Docker%20Desktop%20images_02.png?raw=true)

    ![](https://github.com/chingyaousf/Introduction-to-Building-Docker-for-BERLIN-Pipeline/blob/main/data/BATCH-FLEX-Shiny%20Docker%20Desktop%20images_03.png?raw=true)

-   Run the image using Windows cmd

    command:

    ```         
    docker run -it -p 3838:3838 batch_flex_shinyapp:latest
    ```

Once the container is running, you should be able to access your Shiny app by navigating to:

 - <http://localhost:3838>

 - <http://127.0.0.1:3838>

in your web browser on the host machine.

**4. Use Docker Compose for Building and Running**:

-   Create a Docker Compose YAML file in the Docker project directory (already included in the Docker file).

-   Inside the YAML file, specify the services, image names, ports, and other configurations.

-   Use **`docker-compose build`** to build the image based on the YAML file. You can run this command on Windows command prompt whenever the Dockerfile is modified to update the image.

-   Use **`docker-compose up`** to run the container based on the built image. This is also used to start the container whenever needed.

-   To remove the container and associated resources, use **`docker-compose down`**.

The Docker container is configured to run Shiny app, ensuring that you can efficiently explore the intricacies of BATCH-FLEX-Shiny with ease and consistency. This setup allows for convenient access to Shiny app. Remember that you can run `docker-compose build` or `docker-compose up` every time you modify the Dockerfile to update the image and keep your analysis environment up to date. Please refer to the **Docker Desktop Learning Center** and **Docker Compose documentation** for comprehensive guidance on using these tools effectively.
