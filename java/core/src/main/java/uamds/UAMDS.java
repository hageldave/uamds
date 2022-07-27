package uamds;

import static uamds.other.Utils.sq;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import uamds.optimization.generic.numerics.MatCalc;
import uamds.optimization.generic.problem.ScalarFN;
import uamds.optimization.generic.problem.VectorFN;
import uamds.optimization.generic.solver.GradientDescent;
import uamds.optimization.generic.solver.StochasticGradientDescent;
import uamds.other.NRV;
import uamds.other.NRVSet;
import uamds.other.Ref;

/**
 * Class for performing uncertainty-aware multidimensional scaling on sets of
 * normally distributed random vectors (multivariate normal distributions).
 * 
 * @param <M> matrix data type
 */
public class UAMDS<M> {
	
	protected final MatCalc<M> mc;
	protected boolean stochasticGD = false;
	public boolean verbose = false;
	public final GradientDescent<M> gd;
	public final int lowDim;
	
	/**
	 * Creates a new UAMDS instance that is using the specified matrix calculator object.
	 *  
	 * @param mc matrix calculator
	 * @param lowDim number of dimensions in projection space (usually 2)
	 */
	public UAMDS(MatCalc<M> mc, int lowDim) {
		this.mc = mc;
		this.lowDim = lowDim;
		this.gd = new StochasticGradientDescent<>(mc);
		this.gd.terminationStepSize = 1e-13;
		this.gd.lineSearchFactor = 1e-3;
	}
	
	public UAMDS(MatCalc<M> mc) {
		this(mc,2);
	}
	
	/**
	 * @return true when stochastic gradient descend is used.
	 */
	public boolean isStochasticGDEnabled() {
		return stochasticGD;
	}
	
	/**
	 * En-/disables the use of stochastic gradient descent (SGD).
	 * In each step of SGD, one distribution of the data set is randomly selected, 
	 * and only stress and gradient components that correspond to it are evaluated.
	 * This makes descent steps cheaper to compute (O(n) instead O(n*n)), but more steps must be taken so that all
	 * parts of the stress are considered eventually. 
	 * Overall using SGD usually provides a considerable speed-up.
	 * 
	 * @param stochasticGD true when enabling (disabled by default)
	 */
	public void setStochasticGDEnabled(boolean stochasticGD) {
		this.stochasticGD = stochasticGD;
	}
	
	public NRVSet<M> calculateProjection(NRVSet<M> data, M[][] init, Ref<M[][]> result) {
		return calculateProjection(data, init, result, 100);
	}
	
	public NRVSet<M> calculateProjection(NRVSet<M> data, M[][] init, Ref<M[][]> result, int numDescentSteps) {
		return calculateProjection(data, init, result, numDescentSteps, null);
	}
	
	public NRVSet<M> calculateProjection(NRVSet<M> data, M[][] init, Ref<M[][]> result, int numDescentSteps, Ref<double[][]> loss){
		return calculateProjection(data, init, result, numDescentSteps, loss, null);
	}
	
	/**
	 * Performs a number of iterations of UAMDS.
	 * For each high-dimensional normal distribution N(μ_i, Σ_i) in the data set, its low dimensional
	 * projection is calculated through an affine transform which UAMDS optimizes for.
	 * <p><pre>
	 * N(μ_i, Σ_i) → projection → N(c_i, W_i)
	 * c_i = P_i * μ_i + t_i   (with projection matrix P_i and translation vector t_i)
	 * W_i = P_i * Σ_i * P_iT  (T for transpose)
	 * </pre><p>
	 * However, the algorithm does not compute P_i and t_i, but instead computes
	 * c_i directly and a projection matrix B_i that relates to P_i like so:
	 * <p><pre>
	 * Σ_i = U_i * S_i * V_iT = U_i * S_i * U_iT  (singular value decomposition of covariance matrix)
	 * B_i = P_i * U_i  so that  W_i = B_i * S_i * B_iT  and  P_i = B_i * U_iT
	 * </pre>
	 * 
	 * @param data set of NRVs (normal distributions) N(μ_i, Σ_i)
	 * 
	 * @param init (optional, can be null) initial guess of the affine transforms B_i,c_i for each NRV of the data set.<br>
	 * init[0][i]=B[i], init[1][i]=c[i].<br>
	 * When init=null, random values will be used to initialize.
	 * 
	 * @param result (optional, can be null) will hold the optimized affine transfroms 
	 * B_i,c_i and related P_i,t_i when the method returns.<br>
	 * result.r[0][i]=B[i], result.r[1][i]=c[i], result.r[2][i]=P[i], result.r[3][i]=t[i].
	 * 
	 * @param numDescentSteps maximum number of steps gradient descend will take before returning.
	 * 
	 * @param loss (optional, can be null) will hold the pairwise stress between projected NRVs when the method returns.<br>
	 * loss[i][j]=loss[j][i] = 'stress between distribution i and j in the projection'
	 * 
	 * @param stressComps (optional, can be null) will hold the components 1,2,3 of the pairwise stress between projected NRVs when the method returns.<br>
	 * stress[k][i][j]=loss[k][j][i] = 'stress component k (in 0,1,2) between distribution i and j in the projection'
	 * 
	 * @return projected NRVs (normal distributions) N(c_i, W_i)
	 */
	public NRVSet<M> calculateProjection(NRVSet<M> data, M[][] init, Ref<M[][]> result, int numDescentSteps, Ref<double[][]> loss, Ref<double[][][]> stressComps) {
		int hiDim = data.get(0).d;
		int loDim = lowDim;
		PreCalculatedValues<M> pre = new PreCalculatedValues<>(mc,data);
		
		M[][] solution = 
				stochasticGD 
				? 
				optimizeUAMDS_stochastic(loDim, pre, init, numDescentSteps) 
				: 
				optimizeUAMDS(loDim, pre, init, numDescentSteps);
		// solution extraction and projection
		M[] B = solution[0];
		M[] c = solution[1];
		M[] P = mc.matArray(data.size());
		M[] t = mc.matArray(data.size());
		NRVSet<M> lowPointset = new NRVSet<>();
		for(int i=0; i<data.size(); i++) {
			NRV<M> projected = new NRV<M>(mc, c[i], mc.mult_abcT(B[i], pre.S[i], B[i]));
			lowPointset.add(projected);
			P[i] = mc.mult_abT(B[i], pre.U[i]);
			t[i] = mc.sub(c[i], mc.mult_ab(P[i],pre.mu[i]));
		}
		if(result != null) {
			M[][] resultTransforms = mc.matArray(4, 0);
			resultTransforms[0] = B;
			resultTransforms[1] = c;
			resultTransforms[2] = P;
			resultTransforms[3] = t;
			result.set(resultTransforms);
		}
		
		double[][] loss_ij = new double[data.size()][data.size()];
		double[][][] loss_kij = new double[3][data.size()][data.size()];
		double stress = stressFromProjecton(pre, B, c, hiDim, loDim, loss_ij, loss_kij);
		if(loss!=null)
			loss.set(loss_ij);
		if(stressComps!=null) {
			stressComps.set(loss_kij);
		}
		if(verbose)
			System.out.println("stress=" + stress);
	
		return lowPointset;
	}
	
	protected M[][] optimizeUAMDS(final int loDim, PreCalculatedValues<M> pre, M[][] init, int numDescentSteps) {		
		final int n = pre.n; // number of distributions
		final int hiDim = mc.numElem(pre.mu[0]);
		
		/* initialize affine transforms */
		M[] B,c;
		if(init != null && init.length >= 2 && init[0].length == n) {
			B = Arrays.stream(init[0]).map(mc::copy).toArray(mc::matArray);
			c = Arrays.stream(init[1]).map(mc::copy).toArray(mc::matArray);
			// check
			for(int i=0;i<B.length; i++) {
				if(mc.numCols(B[i])!= hiDim || mc.numRows(B[i])!= loDim || mc.numElem(c[i])!= loDim)
					throw new IllegalArgumentException("something is not matching up with dimensions of provided init");
			}
		} else {
			B = mc.matArray(n);
			c = mc.matArray(n);
			for(int i=0;i<B.length; i++) {
				B[i] = mc.sub(mc.rand(loDim, hiDim),0.5);
				c[i] = mc.sub(mc.rand(loDim),0.5);
			}
		}
		
		/* create the optimization problem (loss function for stress minimization).
		 * This needs objects x, f(x), f'(x)
		 */
		M x = vectorizeAffineTransforms(B, c, null);
		ScalarFN<M> fx = new ScalarFN<M>() {
			@Override
			public double evaluate(M vec) {
				extractAffineTransforms(B, c, vec);
				return stressFromProjecton(pre, B, c, hiDim, loDim, null, null);
			}
		};
		VectorFN<M> dfxa = new VectorFN<M>() {
			M[] B = mc.matArray(n);
			M[] c = mc.matArray(n);
			M grad = mc.copy(x);
			{ // default constructor
				for(int i=0; i<B.length; i++) {
					B[i] = mc.zeros(loDim, hiDim);
					c[i] = mc.zeros(loDim);
				}
			}
			
			@Override
			public M evaluate(M vec) {
				extractAffineTransforms(B, c, vec);
				M[][] dc_dB = gradientFromProjection(pre, B, c, hiDim, loDim);
				return vectorizeAffineTransforms(dc_dB[1], dc_dB[0], grad);
			}
		};
		
		/* (only for debugging)
		NumericGradient<M> dfxn = new NumericGradient<>(mc, fx);
		dfxn.h = 1e-9;
		M[] tempBs = Arrays.stream(projectionsDistrSpace).map(mc::copy).toArray(mc::matArray);
		M[] tempcs = Arrays.stream(loMeans).map(mc::copy).toArray(mc::matArray);
		VectorFN<M> dfx = (x_)->{ // analytic gradient with on the fly numeric gradient check
			var dxCheck = dfxn.central.evaluate(x_);
			var dx_ = dfxa.evaluate(x_);
			var diffGrad = mc.sub(dx_, dxCheck);
			copyFromVectorizedAffineTransform(tempBs, tempcs, diffGrad);
			var db = tempBs;
			var dc = tempcs;
			double diff = mc.norm(diffGrad);
			double lenratio = mc.norm(dx_)/mc.norm(dxCheck);
			double dot = mc.inner(dx_,dxCheck)/(mc.norm(dx_)*mc.norm(dxCheck));
			if(diff > 0.01)
				System.err.println("gradient mismatch! diff="+diff+" dot="+dot+ " lenRatio="+lenratio);
			return dx_;
		};
		*/
		
		// minimizing with gradient descent
		gd.maxDescentSteps = numDescentSteps;
		M xMin = gd.arg_min(fx, dfxa, x, null);
		if(verbose)
			System.out.println("stepsize on termination:"+gd.stepSizeOnTermination);
		
		// solution extraction
		extractAffineTransforms(B, c, xMin);		
		M[][] result = mc.matArray(2, 0);
		result[0] = B;
		result[1] = c;
		return result;
	}
	
	protected M[][] optimizeUAMDS_stochastic(final int loDim, PreCalculatedValues<M> pre, M[][] init, int numDescentSteps) {
		final int n = pre.n; // number of distributions
		final int hiDim = mc.numElem(pre.mu[0]);
		
		/* initialize affine transforms */
		M[] B,c;
		if(init != null && init.length >= 2 && init[0].length == n) {
			B = Arrays.stream(init[0]).map(mc::copy).toArray(mc::matArray);
			c = Arrays.stream(init[1]).map(mc::copy).toArray(mc::matArray);
			// check
			for(int i=0;i<B.length; i++) {
				if(mc.numCols(B[i])!= hiDim || mc.numRows(B[i])!= loDim || mc.numElem(c[i])!= loDim)
					throw new IllegalArgumentException("something is not matching up with dimensions of provided init");
			}
		} else {
			B = mc.matArray(n);
			c = mc.matArray(n);
			for(int i=0;i<B.length; i++) {
				B[i] = mc.sub(mc.rand(loDim, hiDim),0.5);
				c[i] = mc.sub(mc.rand(loDim),0.5);
			}
		}
		
		/* need to base stress and gradient on a random number which decides on the evaluated member i */
		Ref<Integer> currentRandom = new Ref<>(0);
		/* create the optimization problem (loss function for stress minimization).
		 * This needs objects x, f(x), f'(x)
		 */
		M x = vectorizeAffineTransforms(B, c, null);
		ScalarFN<M> fx = new ScalarFN<M>() {
			@Override
			public double evaluate(M vec) {
				extractAffineTransforms(B, c, vec);
				int i=currentRandom.get()%n;
				return stressFromProjection_i(i, pre, B, c, hiDim, loDim, null, null);
			}
		};
		VectorFN<M> dfxa = new VectorFN<M>() {
			M[] B = mc.matArray(n);
			M[] c = mc.matArray(n);
			M grad = mc.copy(x);
			
			M cZero;
			M BZero;
			{ // default constructor
				for(int i=0; i<B.length; i++) {
					B[i] = mc.zeros(loDim, hiDim);
					c[i] = mc.zeros(loDim);
				}
				cZero = mc.zeros(loDim);
				BZero = mc.zeros(loDim, hiDim);
			}
			
			M[] ntimes(int n, M m) {
				M[] arr = mc.matArray(n);
				Arrays.fill(arr, m);
				return arr;
			}
			
			@Override
			public M evaluate(M vec) {
				extractAffineTransforms(B, c, vec);
				
				int i=currentRandom.get()%n;
				M[][] dc_dB_i = gradientFromProjection_i(i, pre, B, c, hiDim, loDim);
				M[] dc = ntimes(n, cZero), dB = ntimes(n, BZero);
				dc[i] = dc_dB_i[0][0];
				dB[i] = dc_dB_i[1][0];
				
				return vectorizeAffineTransforms(dB, dc, grad);
			}
		};
		
		// minimizing with gradient descent
		gd.maxDescentSteps = numDescentSteps*n;
		((StochasticGradientDescent<M>)gd).randRef = currentRandom;
		M xMin = gd.arg_min(fx, dfxa, x, null);
		if(verbose)
			System.out.println("stepsize on termination:"+gd.stepSizeOnTermination);
		
		// solution extraction
		extractAffineTransforms(B, c, xMin);		
		M[][] result = mc.matArray(2, 0);
		result[0] = B;
		result[1] = c;
		return result;
	}
	
	
	protected M vectorizeAffineTransforms(M[] Bs, M[] cs, M target) {
		int hiDim = mc.numCols(Bs[0]);
		int loDim = mc.numRows(Bs[0]);
		int n = Bs.length;
		int nComps = n*hiDim*loDim + n*loDim;
		if(target == null)
			target = mc.zeros(nComps);
		if(mc.numElem(target) != nComps)
			throw new IllegalArgumentException("provided target vector is not of correct size. Should be " + nComps + " but is " + mc.numElem(target));
		
		int i=0;
		for(int j=0; j<Bs.length; j++) {
			M B = Bs[j];
			M c = cs[j];
			mc.copyValues(B, 0, target, i, mc.numElem(B));
			i+= mc.numElem(B);
			mc.copyValues(c, 0, target, i, mc.numElem(c));
			i+= mc.numElem(c);
		}
		return target;
	}
	
	protected void extractAffineTransforms(M[] Bs, M[] cs, M source) {
		int hiDim = mc.numCols(Bs[0]);
		int loDim = mc.numRows(Bs[0]);
		int n = Bs.length;
		int nComps = n*hiDim*loDim + n*loDim;
		if(mc.numElem(source) != nComps)
			throw new IllegalArgumentException("provided target vector is not of correct size. Should be " + nComps + " but is " + mc.numElem(source));
		
		int i=0;
		for(int j=0; j<Bs.length; j++) {
			M B = Bs[j];
			M c = cs[j];
			mc.copyValues(source, i, B, 0, mc.numElem(B));
			i+= mc.numElem(B);
			mc.copyValues(source, i, c, 0, mc.numElem(c));
			i+= mc.numElem(c);
		}
	}


	protected double stressFromProjecton(PreCalculatedValues<M> pre, M[] B, M[] c, int hiDim, int loDim, double[][] loss_ij, double[][][] loss_kij) {
		double sum = 0;
		for(int i=0; i<B.length; i++) {
			for(int j=i; j<B.length; j++) {
				double loss = stressFromProjecton_ij(i, j, pre, B, c, hiDim, loDim, loss_ij, loss_kij);
				sum += loss;
			}
		}
		return sum;
	}
	
	protected double stressFromProjection_i(int i, PreCalculatedValues<M> pre, M[] B, M[] c, int hiDim, int loDim, double[][] loss_ij, double[][][] loss_kij) {
		double sum=0;
		for(int j=0; j<B.length; j++) {
			sum += j<i ? stressFromProjecton_ij(j, i, pre, B, c, hiDim, loDim, loss_ij, loss_kij) : stressFromProjecton_ij(i, j, pre, B, c, hiDim, loDim, loss_ij, loss_kij);
		}
		return sum;
	}
	
	protected double stressFromProjecton_ij(int i, int j, PreCalculatedValues<M> pre, M[] B, M[] c, int hiDim, int loDim, double[][] loss_ij, double[][][] loss_kij) {
		M Si = pre.S[i];
		M Ssqrti = pre.Ssqrt[i];
		M Bi = B[i];
		M ci = c[i];

		M Sj = pre.S[j];
		M Ssqrtj = pre.Ssqrt[j];
		M Bj = B[j];
		M cj = c[j];

		M temp;

		double term1;
		{
			// term 1 part 1 : ||Si - Si^(1/2)Bi^T BiSi^(1/2)||_F^2
			temp = mc.mult_ab(Bi, Ssqrti);
			temp = mc.sub(Si, mc.mult_aTb(temp, temp));
			double part1 = mc.frob2(temp);

			// term 1 part 2 : same as part 1 but with j
			temp = mc.mult_ab(Bj, Ssqrtj);
			temp = mc.sub(Sj, mc.mult_aTb(temp, temp));
			double part2 = mc.frob2(temp);

			// term 1 part 3 :
			var tempStatic = pre.SsqrtiUiTUjSsqrtj[i][j];// mc.mult_aTb(mc.mult_ab(Ui, Ssqrti), mc.mult_ab(Uj, Ssqrtj));
			temp = mc.mult_aTb(mc.mult_ab(Bi, Ssqrti), mc.mult_ab(Bj, Ssqrtj));
			temp = mc.sub(tempStatic, temp);
			double part3 = mc.frob2(temp);

			// term 1
			term1 = 2*(part1+part2)+4*part3;
		}

		double term2;
		{
			// term2 part 1 : sum_k^n [ Si_k * ( <Ui_k, mui-muj> - <Bi_k, ci-cj> )^2 ]
			temp = pre.muisubmujTUi[i][j]; // mc.mult_aTb(mc.sub(mui, muj), Ui);
			temp = mc.sub(temp, mc.mult_aTb(mc.sub(ci, cj), Bi));
			temp = mc.elmmul(temp, temp);
			double part1 = mc.sum(mc.mult_ab(temp, Si));

			temp = pre.muisubmujTUj[i][j]; // mc.mult_aTb(mc.sub(mui, muj), Uj);
			temp = mc.sub(temp, mc.mult_aTb(mc.sub(ci, cj), Bj));
			temp = mc.elmmul(temp, temp);
			double part2 = mc.sum(mc.mult_ab(temp, Sj));

			// term 2
			term2 = part1 + part2;
		}

		double term3;
		{
			double norm1 = pre.norm2muisubmuj[i][j]; // mc.frob2(mc.sub(mui, muj));
			double norm2 = mc.norm2(mc.sub(ci, cj));
			double part1 = norm1-norm2;

			double part2=0, part3=0;
			for(int k=0; k<hiDim; k++) {
				double sigmai = mc.get(Si,k,k);
				M bik = mc.getCol(Bi, k);
				double sigmaj = mc.get(Sj, k,k);
				M bjk = mc.getCol(Bj, k);
				part2 += (1.0 - mc.norm2(bik))*sigmai;
				part3 += (1.0 - mc.norm2(bjk))*sigmaj;
			}

			term3 = sq(part1 + part2 + part3);
		}
		double loss = term1+term2+term3;
		if(loss_ij != null) {
			loss_ij[i][j] = loss_ij[j][i] = loss;
		}
		if(loss_kij != null) {
			loss_kij[0][i][j] = loss_kij[0][j][i] = term1;
			loss_kij[1][i][j] = loss_kij[1][j][i] = term2;
			loss_kij[2][i][j] = loss_kij[2][j][i] = term3;
		}
		return loss;
	}
	
	protected M[][] gradientFromProjection(PreCalculatedValues<M> pre, M[] B, M[] c, int hiDim, int loDim){
		// variables for derivatives w.r.t. c and B
		M[] dc = mc.matArray(c.length);
		M[] dB = mc.matArray(B.length);
		for(int i=0; i<dc.length; i++) {
			dc[i] = mc.zeros(loDim);
			dB[i] = mc.zeros(loDim, hiDim);
		}
		
		for(int i=0; i<B.length; i++) {
			for(int j=i; j<B.length; j++) {
				M[][] g = gradientFromProjection_ij(i, j, pre, B, c, hiDim, loDim);
				mc.add_inp(dc[i], g[0][0]);
				mc.add_inp(dc[j], g[0][1]);
				mc.add_inp(dB[i], g[1][0]);
				mc.add_inp(dB[j], g[1][1]);
			}
		}
		
		M[][] result = mc.matArray(2,0);
		result[0]=dc; result[1]=dB;
		return result;
	}
	
	protected M[][] gradientFromProjection_i(int i, PreCalculatedValues<M> pre, M[] B, M[] c, int hiDim, int loDim){
		M dBi = mc.zeros(loDim, hiDim);
		M dci = mc.zeros(loDim);
		for(int j=0; j<B.length; j++) {
			if(j < i) {
				M[][] g = gradientFromProjection_ij(j, i, pre, B, c, hiDim, loDim);
				mc.add_inp(dci, g[0][1]);
				mc.add_inp(dBi, g[1][1]);
			} else {
				M[][] g = gradientFromProjection_ij(i, j, pre, B, c, hiDim, loDim);
				mc.add_inp(dci, g[0][0]);
				mc.add_inp(dBi, g[1][0]);
			}
		}
		M[][] result = mc.matArray(2, 1);
		result[0][0] = dci;
		result[1][0] = dBi;
		return result;
	}
	
	protected M[][] gradientFromProjection_ij(int i, int j, PreCalculatedValues<M> pre, M[] B, M[] c, int hiDim, int loDim){	
		M mui = pre.mu[i];
		M Si = pre.S[i];
		M Bi = B[i];
		M BiSi = mc.mult_ab(Bi, Si);
		M ci = c[i];
		M dBi = mc.zeros(loDim, hiDim);
		M dci = mc.zeros(loDim);

		M muj = pre.mu[j];
		M Sj = pre.S[j];
		M Bj = B[j];
		M BjSj = mc.mult_ab(Bj, Sj);
		M cj = c[j];
		M dBj = mc.zeros(loDim, hiDim);
		M dcj = mc.zeros(loDim);

		M cisubcj = mc.sub(ci,cj);
		//TODO: precompute mu diff
		M muisubmuj = mc.sub(mui, muj);


		//double term1;
		{
			M Zij = pre.Zij[i][j];// = mc.mult_aTb(Ui, Uj);

			M part1i = mc.sub(mc.mult_ab(mc.mult_abT(BiSi, Bi), BiSi) , mc.mult_ab(BiSi, Si));
			M part1j = mc.sub(mc.mult_ab(mc.mult_abT(BjSj, Bj), BjSj) , mc.mult_ab(BjSj, Sj));

			M part2i = mc.sub(mc.mult_ab(mc.mult_abT(BjSj, Bj), BiSi), mc.mult_ab(mc.mult_abT(BjSj, Zij), Si));
			M part2j = mc.sub(mc.mult_ab(mc.mult_abT(BiSi, Bi), BjSj), mc.mult_ab(mc.mult_ab (BiSi, Zij), Sj));

			mc.add_inp(dBi, mc.scale(mc.add(part1i, part2i), 8));
			mc.add_inp(dBj, mc.scale(mc.add(part1j, part2j), 8));
		}

		//double term2;
		if(i != j){
			M part3i = 
					mc.mult_ab(
							mc.sub( 
									mc.mult_ab(cisubcj, mc.mult_aTb(cisubcj,Bi)) ,  
									mc.mult_ab(cisubcj, pre.muisubmujTUi[i][j]/*mc.mult_aTb(muisubmuj,Ui)*/) ),
							Si);
			M part3j = 
					mc.mult_ab(
							mc.sub( 
									mc.mult_ab(cisubcj, mc.mult_aTb(cisubcj,Bj)) ,  
									mc.mult_ab(cisubcj, pre.muisubmujTUj[i][j]/*mc.mult_aTb(muisubmuj,Uj)*/) ),
							Sj);

			mc.add_inp(dBi, mc.scale(part3i, 2));
			mc.add_inp(dBj, mc.scale(part3j, 2));

			M part4i = 
					mc.mult_ab(
							mc.mult_ab(Bi,Si),
							mc.sub(
									pre.muisubmujTUi_T[i][j]/*mc.mult_aTb(Ui, muisubmuj)*/,
									mc.mult_aTb(Bi, cisubcj))
							);
			M part4j = 
					mc.mult_ab(
							mc.mult_ab(Bj,Sj),
							mc.sub(
									pre.muisubmujTUj_T[i][j]/*mc.mult_aTb(Uj, muisubmuj)*/,
									mc.mult_aTb(Bj, cisubcj))
							); 
			M part4 = mc.scale(mc.add(part4i, part4j), -2);

			mc.add_inp(dci, part4);
			mc.sub_inp(dcj, part4);
		}

		double term3;
		{
			double norm1 = mc.norm2(muisubmuj);
			double norm2 = mc.norm2(cisubcj);
			double part1 = norm1-norm2;

			double part2=0, part3=0;
			for(int k=0; k<hiDim; k++) {
				double sigmai = mc.get(Si,k,k);
				M bik = mc.getCol(Bi,k);
				double sigmaj = mc.get(Sj,k,k);
				M bjk = mc.getCol(Bj,k);
				part2 += (1-mc.norm2(bik))*sigmai;
				part3 += (1-mc.norm2(bjk))*sigmaj;
			}

			term3 = -4*(part1 + part2 + part3);

			mc.add_inp(dBi, mc.scale(BiSi,term3));
			mc.add_inp(dBj, mc.scale(BjSj,term3));

			if(i!=j) {
				mc.add_inp(dci, mc.scale(cisubcj, term3));
				mc.sub_inp(dcj, mc.scale(cisubcj, term3));
			}
		}
		
		M[][] result = mc.matArray(2,2);
		result[0][0]=dci;
		result[1][0]=dBi;
		result[0][1]=dcj;
		result[1][1]=dBj;
		return result;
	}
	
	/**
	 * Class storing constant matrix matrix/vector products occurring
	 * in the stress and its gradient.
	 * 
	 * @param <M> matrix data type
	 */
	public static class PreCalculatedValues<M> {
		public final int n;
		public final M[] U;
		public final M[] S;
		public final M[] Ssqrt;
		public final M[] mu;
		
		public final M[][] SsqrtiUiTUjSsqrtj;
		public final M[][] muisubmujTUi;
		public final M[][] muisubmujTUj;
		public final M[][] muisubmujTUi_T;
		public final M[][] muisubmujTUj_T;
		public final double[][] norm2muisubmuj;
		public final M[][] Zij;


		public PreCalculatedValues(MatCalc<M> mc, NRVSet<M> data) {
			this.n=data.size();
			this.mu = data.stream().map(nrv->nrv.mean).toArray(mc::matArray);
			/* singular value decompositions of covariances */
			List<M[]> svds = data.stream().map(nrv->mc.svd(nrv.cov, true)).collect(Collectors.toList());
			this.U = mc.matArray(n);
			this.S = mc.matArray(n);
			this.Ssqrt = mc.matArray(n);
			for(int i=0; i<n; i++) {
				U[i] = svds.get(i)[0];
				S[i] = svds.get(i)[1];
				Ssqrt[i] = mc.sqrt_inp(mc.copy(S[i]));
			}
			/* storage for precomputed constant matrices */
			this.SsqrtiUiTUjSsqrtj = mc.matArray(n,n);
			this.muisubmujTUi = mc.matArray(n,n);
			this.muisubmujTUj = mc.matArray(n,n);
			this.muisubmujTUi_T = mc.matArray(n,n);
			this.muisubmujTUj_T = mc.matArray(n,n);
			this.norm2muisubmuj = new double[n][n];
			this.Zij = mc.matArray(n,n);
			
			for(int i=0; i<n; i++) {
				M mui = mu[i];
				M Ui = U[i];
				M Ssqrti = Ssqrt[i];
				for(int j=i; j<n; j++) {
					M muj = mu[j];
					M Uj = U[j];
					M Ssqrtj = Ssqrt[j];
					
					{ // precalculation
						SsqrtiUiTUjSsqrtj[i][j] = mc.mult_aTb(mc.mult_ab(Ui, Ssqrti), mc.mult_ab(Uj, Ssqrtj));
						muisubmujTUi[i][j] = mc.mult_aTb(mc.sub(mui, muj), Ui);
						muisubmujTUj[i][j] = mc.mult_aTb(mc.sub(mui, muj), Uj);
						muisubmujTUi_T[i][j] = mc.trp(muisubmujTUi[i][j]);
						muisubmujTUj_T[i][j] = mc.trp(muisubmujTUj[i][j]);
						norm2muisubmuj[i][j] = mc.norm2(mc.sub(mui, muj));
						Zij[i][j] = mc.mult_aTb(Ui, Uj);
					}
				}
			}
		}
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

}
