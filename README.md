# UAMDS
Uncertainty-Aware Multidimensional Scaling (UAMDS) is a dimensionality reduction method for uncertain data. Paper publication available at https://doi.org/10.1109/TVCG.2022.3209420 (to be presented at IEEE VIS 2022).

To model the uncertainty of each element in a dataset, an element is no longer expressed as a vector, but as a random vector with its own multivariate probability distribution.
So instead of reducing a set of high-dimensional data points to low-dimensional points $p_i \rightarrow x_i$, it reduces high-dimensional random vectors to low-dimensional random vectors $P_i \rightarrow X_i$.  
This implementation supports datasets of normally distributed random vectors, which means that each element is a tuple consisting of a mean vector and covariance matrix $P_i \sim N(\mu_i, \Sigma_i)$.

UAMDS computes an optimal affine transformation (linear projection & translation) $\Phi_i(p) = M_i~p + d_i$ for each random vector in the dataset, that transforms the random vector to low-dimensional space, $X_i = \Phi_i(P_i)$. Thus, each low-dimensional random vector will also be normally distributed, $X_i \sim N(M_i ~ \mu_i + d_i, ~ M_i ~ \Sigma_i ~ M_i^\top )$ with projection matrix $M_i$ and translation $d_i$.


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


