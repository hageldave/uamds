# UAMDS
Uncertainty-Aware Multidimensional Scaling

## Java Implementation
The java project consists of a library which is served as a [maven](https://maven.apache.org/what-is-maven.html) artifact.

### Compile & Run Demo
To compile the library you'll need a [Java development kit](https://adoptopenjdk.net/) (version 11 or higher) and either an installation of maven, or an up to date Java IDE with integrated maven support (such as [Eclipse](https://www.eclipse.org/)).
In Ubuntu Linux the former option is quite straight forward.
```sh
sudo apt install openjdk-11-jdk
sudo apt install maven
# cd into java/core and then execute maven build&run:
mvn exec:java -D"exec.mainClass"="demo.Example"
```


