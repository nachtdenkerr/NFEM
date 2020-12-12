package plasticity;

import java.util.ArrayList;
import java.util.List;

import fem_1D.GaussPoint;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVectorRO;
import iceb.plotutils.DynamicDataset;
import iceb.plotutils.ui.PlotPanel;
import linalg.Array1DVector;
import linalg.Array2DMatrix;

public class VMLH2D extends PlasticModel{
	private GaussPoint gp;
	private double sigma_y0, H, K, G;
	private IVectorRO SOID;
	private IMatrixRO II,DEVPRJ,Ce;
		
	private List<Double> eps_plot_list = new ArrayList<Double>();
	private List<Double> q_list = new ArrayList<Double>();
	private List<Boolean> LALGVA = new ArrayList<Boolean>();
	private double DGAMA;
	private static double TOL = 1e-6;
	
	public VMLH2D (GaussPoint gp, double sigma_y, double H, double E, double nuy){
		this.gp = gp; this.H = H; 
		this.sigma_y0 = sigma_y;

		// initialize first history variables
		LALGVA.add(false);
		gp.getMaterial().setSigma(new Array1DVector(0.0, 0.0, 0.0, 0.0));
		
		eps_plot_list.add(0.0);
		q_list.add(0.0);
		
		K = E/(3*(1 - 2*nuy));
		G = E/(2*(1 + nuy));
		
		SOID = new Array1DVector(1.0, 1.0, 0.0, 1.0);

		II = new Array2DMatrix(4,4);
		((Array2DMatrix) II).set(0, 0, 1.0);
		((Array2DMatrix) II).set(1, 1, 1.0);
		((Array2DMatrix) II).set(2, 2, 0.5);
		((Array2DMatrix) II).set(3, 3, 1.0);
		
		DEVPRJ = new Array2DMatrix(4,4);
		for (int a = 0; a <4; a++){
			for (int b = 0; b <4; b++){
				((Array2DMatrix) DEVPRJ).set(a,b,II.get(a, b) - (1.0/3.0)*SOID.get(a)*SOID.get(b));
			}
		}
		Ce = new Array2DMatrix(4,4);
		for (int m = 0; m < 4; m++){
			for (int n = 0; n < 4; n++){
				((Array2DMatrix) Ce).set(m,n,DEVPRJ.get(m, n)*2.0*G 
						+ K*SOID.get(m)*SOID.get(n));
			}
		}
		for (int n = 0; n < 3; n++){
			for (int m = n+1; m < 4; m++){
				((Array2DMatrix) Ce).set(m,n,Ce.get(n, m));
			}
		}
	}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public List<IVectorRO> returnMapping(IVectorRO EPSN1, IVectorRO DEPS, IVectorRO EPSEN, double EPBARN, IVectorRO STRESN) {
		
		double q_plot = 0.0;
		// 1. Get history variables
		boolean IFPLAS = LALGVA.get(0);
		DGAMA = 0.0;
		double SIGMAY = gp.getMaterial().PLFUN(sigma_y0, H, EPBARN);

		List result = new ArrayList();
		result.add(EPSN1); 
		
		// 2. Compute elastic trial-strains
		IVectorRO STRAT = EPSEN.add(DEPS);
//		System.out.println("STRAT:" + MatrixFormat.format(STRAT));
		
		// Volumetric strain and pressure stress
		double EEV = STRAT.get(0) + STRAT.get(1) + STRAT.get(3);
		double P = EEV*K; 

		// Elastic trial deviatoric strain
		double EEVD3 = EEV/3.0;
		IVectorRO EET = STRAT.subtract(SOID.multiply(EEVD3));
		((Array1DVector) EET).set(2, STRAT.get(2)/2.0);
//		System.out.println("DEPS:" + MatrixFormat.format(DEPS));

		// 3. Compute trial effective stress and uni-axial yield stress
		double VARJ2T = 0.0;
		VARJ2T = 4.0*G*G*(Math.pow(EET.get(2), 2.0) +
				0.5*(Math.pow(EET.get(0), 2.0) + Math.pow(EET.get(1), 2.0) + (Math.pow(EET.get(3), 2.0)))); 
		
		double QTRIAL = Math.sqrt(3.0*VARJ2T);
		
		// 4. Check for plastic admissibility
		double PHI = QTRIAL - SIGMAY;
		IVectorRO STRES;
		IVectorRO RSTAVA = EPSEN;
		if (PHI / SIGMAY > TOL) {
			// Plastic step: Apply return mapping - use Newton-Raphson algorithm
			// to solve the return mapping equation (Box 7.4)
			
//			System.out.println("!!!Plastic step at Gauss Point " + "!!!");

//			System.out.println("Phi: " + PHI);
			
			IFPLAS = true;
			double EPBAR = EPBARN;
			for (int n = 0; n < 10; n++){
				double DENOM = -3.0*G - gp.getMaterial().DPLFUN(sigma_y0, H, EPBAR);;
				double DDGAMA = -PHI/DENOM;
				DGAMA = DGAMA + DDGAMA;

				EPBAR = EPBAR + DDGAMA;
				SIGMAY = gp.getMaterial().PLFUN(sigma_y0, H, EPBAR);
				PHI = QTRIAL - 3.0*G*DGAMA - SIGMAY;
//				System.out.println("DGAMA iter:" + DGAMA);

				double RESNOR = Math.abs(PHI / SIGMAY);
				if (RESNOR <= TOL){
					
					// update stress components
					double FACTOR = 2.0*G * (1.0 - 3.0*G*DGAMA/QTRIAL);
					STRES = EET.multiply(FACTOR).add(SOID.multiply(P));
					gp.getMaterial().setSigma(STRES);
					result.add(STRES);	
					//System.out.println("EET.multiply(FACTOR):" + MatrixFormat.format(EET.multiply(FACTOR)));

					// compute converged elastic (engineering) strain components
					FACTOR = FACTOR / (2.0*G);
					RSTAVA = EET.multiply(FACTOR).add(SOID.multiply(EEVD3));
					((Array1DVector) RSTAVA).set(2,EET.multiply(FACTOR).get(2)*2.0);
					
//					System.out.println("QTRIAL:" + (QTRIAL));
					
					// update accumulated plastic strain
					result.add(EPBAR);
					gp.getMaterial().setEpsPlastic(EPBAR);

					break;
				}
			}
			q_plot = QTRIAL - 3.0*G*DGAMA;
		}
		else {
			// Elastic step: Update stress using linear elastic law
			STRES = EET.multiply(2.0*G).add(SOID.multiply(P));

			// elastic engineering strain
			RSTAVA = STRAT;
			gp.getMaterial().setSigma(STRES);
			//gp.setDGAMA(DGAMA);
			double EPBAR = EPBARN;
			gp.getMaterial().setEpsPlastic(EPBAR);
			q_plot = QTRIAL;
			result.add(STRES);			
			result.add(EPBAR);
		}

		LALGVA.set(0, IFPLAS);
		q_list.add(q_plot);
		eps_plot_list.add(EPSN1.get(2));
		System.out.println();
		
		result.add(RSTAVA);
		gp.getMaterial().setEpsElastic(RSTAVA);
		return result;
	}
	
	public void tangentUpdate (double EPBAR) {
		
		IVectorRO STRES = gp.getMaterial().getSigma();
		//double DGAMA = gp.getDGAMA();
		IMatrixRO DMATX = new Array2DMatrix(4,4);
		boolean IFPLAS = LALGVA.get(0);

		if (IFPLAS == true){
			// Compute elastoplastic consistent tangent
			double R003D2 = Math.sqrt(3.0/2.0);
			// Hydrostatic pressure
			double P = 1.0/3.0 * (STRES.get(0) + STRES.get(1) + STRES.get(3));	
			// Deviatoric stress components
			IVectorRO S = STRES.subtract(SOID.multiply(P));				
			
			// Recover last elastic trial von Mises effective stress
			double SNORM = Math.pow(S.get(0),2) + Math.pow(S.get(1),2) + Math.pow(S.get(3),2) 
					+ 2.0*Math.pow(S.get(2), 2);
			
			SNORM = Math.sqrt(SNORM);
			double Q = R003D2*SNORM;
			double QTRIAL = Q + 3.0*G*DGAMA;
			
			// Assemble elastoplastic tangent (upper triangle only)
			double AFACT = 2.0*G*(1.0 - 3.0*G*DGAMA/QTRIAL);
			double BFACT = 6.0*G*G*(DGAMA/QTRIAL - 1.0/(3.0*G + gp.getMaterial().DPLFUN(sigma_y0, H, EPBAR)))/(SNORM*SNORM);

			for (int m = 0; m < 4; m++){
				for (int n = m; n < 4; n++){
					((Array2DMatrix) DMATX).set(m,n,DEVPRJ.get(m, n)*AFACT 
							+ S.get(m)*S.get(n)*BFACT + K*SOID.get(m)*SOID.get(n));
				}
			}
		}
		else {
			// Compute elasticity matrix (upper triangle only)
			for (int m = 0; m < 4; m++){
				for (int n = 0; n < 4; n++){
					((Array2DMatrix) DMATX).set(m,n,DEVPRJ.get(m, n)*2.0*G 
							+ K*SOID.get(m)*SOID.get(n));
				}
			}
		}
		// Assemble lower triangle
		for (int n = 0; n < 3; n++){
			for (int m = n+1; m < 4; m++){
				((Array2DMatrix) DMATX).set(m,n,DMATX.get(n, m));
			}
		}
		
		//if (IFPLAS == true){
			//System.out.println("DMATX:\n" + MatrixFormat.format(DMATX));
		//}
		gp.getMaterial().setTangent(DMATX);
	}
	
	public PlotPanel drawPlot(){
		List<Double> e = eps_plot_list;
		List<Double> s = q_list;
		DynamicDataset data = new DynamicDataset();

		for (int i = 0; i < s.size(); i++){
			data.addX(e.get(i));
			data.addY("Load Displacement Curve", s.get(i));
		}
		data.setComplete(true);
		
		PlotPanel plot = new PlotPanel();		
		plot.addDataset(data);
		plot.setXLabel("Strain epsilon");
		plot.setYLabel("Stress sigma");
		//Color c = new Color(220, 220, 220);
		//plot.setBackground(c);

		return plot;
	}
}