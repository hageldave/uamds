# UAMDS
Uncertainty-Aware Multidimensional Scaling (to be presented at IEEE VIS)

## Java Implementation (Reference Implementation)
The java project consists of a library which is served as a [maven](https://maven.apache.org/what-is-maven.html) artifact.

### Compile & Run Demo
To compile the library you'll need a [Java development kit](https://adoptopenjdk.net/) (version 11 or higher) and either an installation of maven, or an up to date Java IDE with integrated maven support (such as [Eclipse](https://www.eclipse.org/)).
In Ubuntu Linux the former option is quite straight forward.
```sh
sudo apt install openjdk-11-jdk
sudo apt install maven
# cd into java/core and then execute maven build&run:
mvn clean compile exec:java -D"exec.mainClass"="uamds.demo.Example"
```

### Use in your own project
The version 0.0.2 (as archived for the paper publication) has to be manually installed. Future releases could be made available through maven central.
For building and installing the maven artifact, the same prerequisites as for running the demo apply. To build and install on CLI run the following:
```sh
# cd into java/core and then execute installation via maven
mvn clean install
```
This will install the artifact into your local maven repository. After that you can use it as a dependency in another maven project.
```xml
<dependency>
	<groupId>de.uni-stuttgart.visus</groupId>
	<artifactId>uamds-core</artifactId>
	<version>0.0.2</version>
</dependency>
```


