package uamds;

import static uamds.Utils.sq;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import optimization.generic.numerics.MatCalc;
import optimization.generic.problem.ScalarFN;
import optimization.generic.problem.VectorFN;
import optimization.generic.solver.GradientDescent;

public class UAMDS<M> {
	
	protected final MatCalc<M> mc;
	public boolean verbose = false;
	
	public UAMDS(MatCalc<M> mc) {
		this.mc = mc;
	}
	
	public RVPointSet<M> calculateProjection(RVPointSet<M> data, M[][] init, Ref<M[][]> result) {
		return calculateProjection(data, init, result, 100);
	}
	
	public RVPointSet<M> calculateProjection(RVPointSet<M> data, M[][] init, Ref<M[][]> result, int numDescentSteps) {
		return calculateProjection(data, init, result, numDescentSteps, null);
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
	 * @return projected NRVs (normal distributions) N(c_i, W_i)
	 */
	public RVPointSet<M> calculateProjection(RVPointSet<M> data, M[][] init, Ref<M[][]> result, int numDescentSteps, Ref<double[][]> loss) {
		int hiDim = data.get(0).d;
		int loDim = 2;
		PreCalculatedValues pre = new PreCalculatedValues(data);
		
		M[][] affineTransforms = optimizeUAMDS(loDim, pre, init, numDescentSteps);
		// solution extraction and projection
		M[] B = affineTransforms[0];
		M[] c = affineTransforms[1];
		M[] projections = mc.matArray(data.size());
		M[] translations = mc.matArray(data.size());
		RVPointSet<M> lowPointset = new RVPointSet<>();
		for(int i=0; i<data.size(); i++) {
			NRV<M> projected = new NRV<M>(mc, c[i], mc.mult_abcT(B[i], pre.S[i], B[i]));
			lowPointset.add(projected);
			projections[i] = mc.mult_abT(B[i], pre.U[i]);
			translations[i] = mc.sub(c[i], mc.mult_ab(projections[i],pre.mu[i]));
		}
		if(result != null) {
			M[][] resultTransforms = mc.matArray(4, 0);
			resultTransforms[0] = B;
			resultTransforms[1] = c;
			resultTransforms[2] = projections;
			resultTransforms[3] = translations;
			result.set(resultTransforms);
		}
		
		double[][] loss_ij = new double[data.size()][data.size()];
		double stress = stressFromProjecton(pre, B, c, hiDim, loDim, loss_ij);
		if(loss!=null)
			loss.set(loss_ij);
		if(verbose)
			System.out.println("stress=" + stress);
	
		return lowPointset;
	}
	
	
	protected M[][] optimizeUAMDS(final int loDim, PreCalculatedValues pre, M[][] init, int numDescentSteps) {		
		final int n = pre.n; // number of distributions
		final int hiDim = mc.numElem(pre.mu[0]);
		
		/* initialize affine transforms */
		M[] projectionsDistrSpace;
		M[] loMeans;
		if(init != null && init.length >= 2 && init[0].length == n) {
			projectionsDistrSpace = Arrays.stream(init[0])
					.map(mc::copy)
					.toArray(mc::matArray);
			loMeans = Arrays.stream(init[1])
					.map(mc::copy)
					.toArray(mc::matArray);
			// check
			for(int i=0;i<projectionsDistrSpace.length; i++) {
				if(mc.numCols(projectionsDistrSpace[i])!= hiDim || mc.numRows(projectionsDistrSpace[i])!= loDim || mc.numElem(loMeans[i])!= loDim)
					throw new IllegalArgumentException("somethings not matching up with dimensions of provided init");
			}
		} else {
			projectionsDistrSpace = mc.matArray(n);
			loMeans = mc.matArray(n);
			for(int i=0;i<projectionsDistrSpace.length; i++) {
				projectionsDistrSpace[i] = mc.sub(mc.rand(loDim, hiDim),0.5);
				loMeans[i] = mc.sub(mc.rand(loDim),0.5);
			}
		}
		
		/* create the optimization problem (loss function for stress minimization).
		 * This needs objects x, f(x), f'(x)
		 */
		M x = vectorizeAffineTransforms(projectionsDistrSpace, loMeans, null);
		ScalarFN<M> fx = new ScalarFN<M>() {
			@Override
			public double evaluate(M vec) {
				extractAffineTransform(projectionsDistrSpace, loMeans, vec);
				return stressFromProjecton(pre, projectionsDistrSpace, loMeans, hiDim, loDim, null);
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
				extractAffineTransform(B, c, vec);
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
		
		// gradient descent
		GradientDescent<M> gd = new GradientDescent<>(mc);
 		gd.maxDescentSteps = numDescentSteps;
		gd.terminationStepSize = 1e-12;
		gd.lineSearchFactor = 1e-3;
		// minimizing
		M xMin = gd.arg_min(fx, dfxa, x, null);
		if(verbose)
			System.out.println("stepsize on termination:"+gd.stepSizeOnTermination);
		
		// solution extraction
		extractAffineTransform(projectionsDistrSpace, loMeans, xMin);		
		M[][] result = mc.matArray(2, 0);
		result[0] = projectionsDistrSpace;
		result[1] = loMeans;
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
	
	protected void extractAffineTransform(M[] Bs, M[] cs, M source) {
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


	protected double stressFromProjecton(PreCalculatedValues pre, M[] B, M[] c, int hiDim, int loDim, double[][] loss_ij) {
		double sum = 0;
		
		for(int i=0; i<B.length; i++) {
			M Si = pre.S[i];
			M Ssqrti = pre.Ssqrt[i];
			M Bi = B[i];
			M ci = c[i];
			for(int j=i; j<B.length; j++) {
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
				sum += loss;
				if(loss_ij != null) {
					loss_ij[i][j] = loss_ij[j][i] = loss;
				}
			}
		}
		
		return sum;
	}
	
	protected M[][] gradientFromProjection(PreCalculatedValues pre, M[] B, M[] c, int hiDim, int loDim){
		// variables for derivatives w.r.t. c and B
		M[] dc = mc.matArray(c.length);
		M[] dB = mc.matArray(B.length);
		for(int i=0; i<dc.length; i++) {
			dc[i] = mc.zeros(loDim);
			dB[i] = mc.zeros(loDim, hiDim);
		}
		// precomputing helper variables
		M[] BS = IntStream.range(0, B.length).mapToObj(i->mc.mult_ab(B[i],pre.S[i])).toArray(mc::matArray);
		
		for(int i=0; i<B.length; i++) {
			M mui = pre.mu[i];
			M Si = pre.S[i];
			M Bi = B[i];
			M BiSi = BS[i];
			M ci = c[i];
			M dBi = dB[i];
			M dci = dc[i];
			
			for(int j=i; j<B.length; j++) {
				M muj = pre.mu[j];
				M Sj = pre.S[j];
				M Bj = B[j];
				M BjSj = BS[j];
				M cj = c[j];
				M dBj = dB[j];
				M dcj = dc[j];
				
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
			}
		}
		
		M[][] result = mc.matArray(2,0);
		result[0]=dc; result[1]=dB;
		return result;
	}
	
	
	private class PreCalculatedValues {
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
		
		
//		public PreCalculatedValues(List<M[]> uss, List<M> mus) {
//			int n=uss.size();
//			this.SsqrtiUiTUjSsqrtj = mc.matArray(n,n);
//			this.muisubmujTUi = mc.matArray(n,n);
//			this.muisubmujTUj = mc.matArray(n,n);
//			this.muisubmujTUi_T = mc.matArray(n,n);
//			this.muisubmujTUj_T = mc.matArray(n,n);
//			this.norm2muisubmuj = new double[n][n];
//			this.Zij = mc.matArray(n,n);
//			
//			for(int i=0; i<n; i++) {
//				M[] uss_i = uss.get(i);
//				M mui = mus.get(i);
//				M Ui = uss_i[0];
//				M Ssqrti = uss_i[2];
//				for(int j=i; j<n; j++) {
//					M[] uss_j = uss.get(j);
//					M muj = mus.get(j);
//					M Uj = uss_j[0];
//					M Ssqrtj = uss_j[2];
//					
//					{ // precalculation
//						SsqrtiUiTUjSsqrtj[i][j] = mc.mult_aTb(mc.mult_ab(Ui, Ssqrti), mc.mult_ab(Uj, Ssqrtj));
//						muisubmujTUi[i][j] = mc.mult_aTb(mc.sub(mui, muj), Ui);
//						muisubmujTUj[i][j] = mc.mult_aTb(mc.sub(mui, muj), Uj);
//						muisubmujTUi_T[i][j] = mc.trp(muisubmujTUi[i][j]);
//						muisubmujTUj_T[i][j] = mc.trp(muisubmujTUj[i][j]);
//						norm2muisubmuj[i][j] = mc.norm2(mc.sub(mui, muj));
//						Zij[i][j] = mc.mult_aTb(Ui, Uj);
//					}
//				}
//			}
//		}


		public PreCalculatedValues(RVPointSet<M> data) {
			this.n=data.size();
			this.mu = data.stream().map(nrv->nrv.mean).toArray(mc::matArray);
			List<M[]> svds = data.stream().map(nrv->mc.svd(nrv.cov, true)).collect(Collectors.toList());
			this.U = mc.matArray(n);
			this.S = mc.matArray(n);
			this.Ssqrt = mc.matArray(n);
			for(int i=0; i<n; i++) {
				U[i] = svds.get(i)[0];
				S[i] = svds.get(i)[1];
				Ssqrt[i] = mc.sqrt_inp(mc.copy(S[i]));
			}
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
