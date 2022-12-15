# UAMDS
Uncertainty-Aware Multidimensional Scaling (UAMDS) is a dimensionality reduction method for uncertain data. The paper publication titled *"Uncertainty-Aware Multidimensional Scaling" by David Hägele, Tim Krake, and Daniel Weiskopf*, is available at https://doi.org/10.1109/TVCG.2022.3209420 (open access, best paper at IEEE VIS 2022).

To model the uncertainty of each element in a dataset, an element is no longer expressed as a vector, but as a random vector with its own multivariate probability distribution.
So instead of reducing a set of high-dimensional data points to low-dimensional points $p_i \rightarrow x_i$, it reduces high-dimensional random vectors to low-dimensional random vectors $P_i \rightarrow X_i$.  
This implementation supports datasets of normally distributed random vectors, which means that each element is a tuple consisting of a mean vector and covariance matrix $P_i \sim N(\mu_i, \Sigma_i)$.

UAMDS computes an optimal affine transformation (linear projection & translation) $\Phi_i(p) = M_i~p + d_i$ for each random vector in the dataset, that transforms the random vector to low-dimensional space, $X_i = \Phi_i(P_i)$. Thus, each low-dimensional random vector will also be normally distributed, $X_i \sim N(M_i ~ \mu_i + d_i, ~ M_i ~ \Sigma_i ~ M_i^\top )$ with projection matrix $M_i$ and translation $d_i$.


## Java Implementation (Reference Implementation)
The java project consists of a library which is served as a [maven](https://maven.apache.org/what-is-maven.html) artifact.   
[![Maven Central](https://img.shields.io/maven-central/v/com.github.hageldave.uamds/uamds-core.svg)](https://search.maven.org/search?q=g:com.github.hageldave.uamds)

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
The maven artifact for version 0.0.3 has been deployed to the central maven repository. So you can simply use it with maven out of the box (maven will download it) using the following dependency.

```xml
<dependency>
	<groupId>com.github.hageldave.uamds</groupId>
	<artifactId>uamds-core</artifactId>
	<version>0.0.3</version>
</dependency>
```
---
The version 0.0.2 (as archived for the paper publication) has to be manually installed. Future releases are made available through maven central.
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

### Replicability Stamp: Instructions to reproduce the teaser figure of the paper from scratch
1. Be on Ubuntu Linux 22.04 (other OS are also possible, but the following instructions are for Ubuntu 22.04)
2. ```sudo apt install openjdk-17-jdk``` installs Java development kit 17 (11 is also fine) 
3. ```sudo apt install maven``` installs Apache Maven
4. ```git clone https://github.com/hageldave/uamds.git``` clones this repository 
5. ```cd uamds/```
6. ```git checkout 15bdd7101bca6f04aae73a98a666e001bf0c184a``` checks out the relevant commit 
7. ```cd java/core/```
8. ```mvn install``` builds and installs the UAMDS core library via Maven
9. ```cd ../plots/``` 
10. ```mvn compile exec:java -D"exec.mainClass"="uamds.plots.Teaser"``` builds and runs code that produces the figure
11. ```CTRL + C``` to terminate the application from command window (alternatively close the Java application window)
12. ```xdg-open teaser.svg``` opens the teaser figure with your default application for svg files, or use ```firefox teaser.svg```

Alternatively you can run the provided scripts like this: ```./install-dependencies.sh; ./replicate-figure.sh```.
The *install-dependencies.sh* script may require super user privileges: ```sudo ./install-dependencies.sh```.
The *replicate-figure.sh* script will clone this repository (again), checkout the required commit, compile the code and run it. A Java application window will open that displays the figure, which will also be exported to *teaser.svg*.


