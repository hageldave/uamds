# UAMDS
Uncertainty-Aware Multidimensional Scaling (UAMDS) is a dimensionality reduction method for uncertain data. The paper publication titled *"Uncertainty-Aware Multidimensional Scaling" by David Hägele, Tim Krake, and Daniel Weiskopf*, is available at https://doi.org/10.1109/TVCG.2022.3209420 (open access, best paper at IEEE VIS 2022).

To model the uncertainty of each element in a dataset, an element is no longer expressed as a vector, but as a random vector with its own multivariate probability distribution.
So instead of reducing a set of high-dimensional data points to low-dimensional points $p_i \rightarrow x_i$, it reduces high-dimensional random vectors to low-dimensional random vectors $P_i \rightarrow X_i$.  
This implementation supports datasets of normally distributed random vectors, which means that each element is a tuple consisting of a mean vector and covariance matrix $P_i \sim N(\mu_i, \Sigma_i)$.

UAMDS computes an optimal affine transformation (linear projection & translation) $\Phi_i(p) = M_i~p + d_i$ for each random vector in the dataset, that transforms the random vector to low-dimensional space, $X_i = \Phi_i(P_i)$. Thus, each low-dimensional random vector will also be normally distributed, $X_i \sim N(M_i ~ \mu_i + d_i, ~ M_i ~ \Sigma_i ~ M_i^\top )$ with projection matrix $M_i$ and translation $d_i$.

---
<picture align="center"><image src="https://github.com/hageldave/uamds/blob/main/images/teaser.svg"/></picture>


## Python Implementation (Port)
A Python port of the algorithm can be found at https://github.com/hageldave/UAMDS-python and as part of https://github.com/hageldave/uadapy. Both projects work but are immature, i.e., lack documentation and are not yet deployed as a package.

## Java Implementation (Reference Implementation)
The project consists of a java library which uses [maven](https://maven.apache.org/what-is-maven.html) as dependency management and build system. 
The library is also served as a maven artifact through the central repository.

[![Maven Central](https://img.shields.io/maven-central/v/com.github.hageldave.uamds/uamds-core.svg)](https://central.sonatype.com/namespace/com.github.hageldave.uamds)

### Compile & Run Demo
To compile the library you'll need a [Java development kit](https://adoptopenjdk.net/) (version 11 or higher) and either an installation of maven, or an up to date Java IDE with integrated maven support (such as [Eclipse](https://www.eclipse.org/)).
In Ubuntu Linux the former option is quite straight forward.
```sh
sudo apt install openjdk-11-jdk
sudo apt install maven
# cd into java/core and then execute maven build&run:
mvn clean compile exec:java -D"exec.mainClass"="uamds.demo.Example"
```

### Integrate in your own project
The maven artifact has been deployed to the [central maven repository](https://repo.maven.apache.org/maven2/com/github/hageldave/uamds/uamds-core/) from version 0.0.3 onwards. So you can simply use it with maven out of the box using the following dependency (maven will download it).

```xml
<dependency>
	<groupId>com.github.hageldave.uamds</groupId>
	<artifactId>uamds-core</artifactId>
	<version>0.1.0</version>
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

### How to use
The implementation uses a generic matrix data type `<M>`. 
This type is determined by the *matrix calculator* object `MatCalc<M>` which wraps around your linear algebra library of choice and offers common operations like creating vectors and matrices, multiplication of matrices, or decompositions like SVD.
This library comes with a wrapper for [EJML's](https://github.com/lessthanoptimal/ejml) `DMatrixRMaj` matrix type. If you want to use another library under the hood, you'll need to write your own implementation of the `MatCalc<M>` interface.
```java
public static void doLinAlg() {
  MatCalc<DMatrixRMaj> mc = new MatCalcEJML();
  DMatrixRMaj zeros = mc.zeros(3);
  DMatrixRMaj ones = mc.add(zeros, 1.0);
  DMatrixRMaj permutation = mc.matOf(new double[][] {{0,1,0},{0,0,1},{1,0,0}});
}
	
public static <M> void doLinAlgGeneric(MatCalc<M> mc) {
  M identity = mc.eye(3);
  M scaling = mc.scale(identity, 1.0/mc.numRows(identity));
  M vec = mc.vecOf(1.0, 2.0, 3.0);
  M result = mc.mult_aTbc(vec, scaling, vec); // a^T * b * c
  double value = mc.get(result, 0, 0);
}
```
To use UAMDS your data set needs to be expressed as a set of normally distributed random vectors. This library has a data type `NRV<M>` for such random vectors and a collection data type `NRVSet<M>`.
A normally distributed vector consists of a mean vector and a covariance matrix.
The `NRV<M>` class offers some common operations like sampling from the distribution or obtaining the corresponding probability density function (PDF).
```java
public static <M> void doDatasetThings(MatCalc<M> mc) {
  int numDims = 3;
  NRV<M> standardNormal = new NRV<>(mc, numDims);
  NRV<M> myNRV = new NRV<>(mc, mc.vecOf(1.0, 2.0, 3.0), mc.eye(numDims, 1.337));
  NRV<M> randCovNRV = new NRV<>(mc, mc.zeros(numDims), NRV.randCov(mc, numDims));
  
  NRVSet<M> dataset = new NRVSet<>();
  dataset.addAll(Arrays.asList(standardNormal, myNRV, randCovNRV));
  M dataSetSamples = dataset.stream()
    // draw 1000 samples form each random vector
    .map(nrv -> nrv.drawSamples(1000))
    // stack all samples to create a large matrix
    .reduce(mc::concatVert)
    .get();
  double[][] samples = mc.toArray2D(dataSetSamples);
}
```
To perform UAMDS an instance of the `UAMDS<M>` class needs to be created. In the simplest case you only want to obtain the projected low dimensional random vectors from your dataset.
The following code will use a random initialization of the affine transformations and perform 100 iterations of gradient descent. 
Then the resulting transformations are used to project the dataset which is then returned.
```java
public static <M> void doUAMDS(MatCalc<M> mc, NRVSet<M> dataset) {
  int lowDims = 2;
  UAMDS<M> uamds = new UAMDS<>(mc, lowDims);
  NRVSet<M> projectedDataset = uamds.calculateProjection(dataset, null, null);
  ...
```
It's likely that 100 iterations are not enough, so you may want to continue iterating. In this case, the resulting affine transformations need to be used as initalization for another invocation of uamds.
To get the affine transformations, a place (reference/pointer `Ref<...>`) needs to be specified where these will be put.
```java
  ...
  Ref<M[][]> affineTransf = new Ref<>();
  NRVSet<M> intermediateResult = uamds.calculateProjection(dataset, null, affineTransf);
  M[][] initialization = affineTransf.get();
  NRVSet<M> refinedResult = uamds.calculateProjection(dataset, initialization, affineTransf);
  ...
```
There are more things that can be tweaked and retrieved, such as the number of iterations to be performed, the resulting stress between pairs of distributions, whether to use stochastic gradient descent or not, hyper parameters for gradient descent, and more.
```java
  ...
  int nIterations = 200;
  Ref<double[][]> stress = new Ref<>();
  Ref<double[][][]> stressDetailed = new Ref<>();
  uamds.setStochasticGDEnabled(true);
  uamds.gd.getHyperparams().set(GradientDescent.PARAM_MAX_LINESEARCH_ITER, 10);
  uamds.calculateProjection(dataset, affineTransf.get(), affineTransf, nIterations, stress, stressDetailed);
}
```
For more clever kinds of initialization the `Initialization` class can be used. 
- There is an improved random intialization that takes the distances within the dataset into account, 
- uncertainty-aware PCA initializer when you want to start from a variance maximized linear projection, 
- and an initialization where UAMDS is performed without variance information, basically mimicking MDS on the means.
```java
M[][] init1 = Initialization.initRandom(mc, dataset);
M[][] init2 = Initialization.initFromUAPCA(mc, dataset);
M[][] init3 = Initialization.initWithoutVariance(mc, dataset);
```



## Replicability Stamp: Instructions to reproduce the teaser figure of the paper from scratch
These are the instructions to replicate the teaser figure (Figure 1) of the paper. Please note that this is the only reproducible figure.
Other figures were created with non-deterministic code (using random intializations) that yields slightly different results on each run and is not part of this repository.

1. Be on Ubuntu Linux 22.04 (other OS are also possible, but the following instructions are for Ubuntu 22.04)
2. ```sudo apt install openjdk-17-jdk``` installs Java development kit 17 (11 is also fine) 
3. ```sudo apt install maven``` installs Apache Maven
4. ```git clone https://github.com/hageldave/uamds.git``` clones this repository 
5. ```cd uamds/```
6. ```git checkout replicability-stamp``` checks out the replicability version (future versions may behave differently)
7. ```cd java/core/```
8. ```mvn install``` builds and installs the UAMDS core library via Maven
9. ```cd ../plots/``` 
10. ```mvn compile exec:java -D"exec.mainClass"="uamds.plots.Teaser"``` builds and runs code that produces the figure
11. ```CTRL + C``` to terminate the application from command window (alternatively close the Java application window)
12. ```xdg-open teaser.svg``` opens the teaser figure with your default application for svg files, or use ```firefox teaser.svg```

Alternatively you can run the provided scripts like this: ```./install-dependencies.sh; ./replicate-figure.sh```.
The *install-dependencies.sh* script may require super user privileges: ```sudo ./install-dependencies.sh```.
The *replicate-figure.sh* script will clone this repository (again), checkout the required commit, compile the code and run it. A Java application window will open that displays the figure, which will also be exported to *teaser.svg*.



