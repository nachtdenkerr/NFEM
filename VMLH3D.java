package plasticity;

//import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import fem_1D.GaussPoint;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVectorRO;
import iceb.plotutils.DynamicDataset;
import iceb.plotutils.ui.PlotPanel;
import linalg.Array1DVector;
import linalg.Array2DMatrix;

public class VMLH3D extends PlasticModel{
	
	private GaussPoint gp;
	private double sigma_y0, H, G, K;
	private IVectorRO SOID;
	private IMatrixRO IxI,II,DEVPRJ;
	
	private List<Double> eps_plot_list = new ArrayList<Double>();
	private List<Double> q_list = new ArrayList<Double>();
	private List<Boolean> LALGVA = new ArrayList<Boolean>();
	private double DGAMA;
	private static double TOL = 1.0e-6;
	
	public VMLH3D (GaussPoint gp, double sigma_y, double H, double E, double nuy){
		this.gp = gp; this.H = H; this.sigma_y0 = sigma_y;
		
		// initialize first history variables
		IVectorRO v0 = new Array1DVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		double c0 = 0.0;

		for (int i = 0; i <8; i++){
			LALGVA.add(false);
		}
		
		gp.getMaterial().setSigma(v0);
		
		
		eps_plot_list.add(c0);
		q_list.add(c0);
		
		K = E/(3*(1 - 2*nuy));
		G = E/(2*(1 + nuy));
		
		IxI = new Array2DMatrix(6,6);
		II = new Array2DMatrix(6,6);

		for (int i = 0; i < 3; i++){
			for (int k = 0; k < 3; k++){
				((Array2DMatrix) IxI).set(i,k,1.0);
			}
			((Array2DMatrix) II).set(i,i,1.0);
			((Array2DMatrix) II).set(i+3,i+3,0.5);
		}
		
		DEVPRJ = II.subtract(IxI.multiply(1.0/3.0));
		SOID = new Array1DVector(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	public List returnMapping(IVectorRO EPSN1, IVectorRO DEPS, IVectorRO EPSEN, double EPBARN, IVectorRO STRESN) {
		double q_plot = 0.0;
		
		List result = new ArrayList();
		result.add(EPSN1); 		
		
		// 1. Get history variables
		boolean IFPLAS = LALGVA.get(0);
		DGAMA = 0.0;
		double SIGMAY = gp.getMaterial().PLFUN(sigma_y0, H, EPBARN);
				
		// 2. Compute elastic trial-strains
		IVectorRO STRAT = EPSEN.add(DEPS);
		
		// Volumetric strain and pressure stress
		double EEV = STRAT.get(0) + STRAT.get(1) + STRAT.get(2);
		double P = EEV*K; 
		
		// Elastic trial deviatoric strain
		double EEVD3 = EEV/3.0;
		IVectorRO EET = STRAT.subtract(SOID.multiply(EEVD3));
		for (int k = 0; k <3; k++){
			((Array1DVector) EET).set(k+3, STRAT.get(k+3)/2.0);
		}
		
		// 3. Compute trial effective stress and uni-axial yield stress
		double VARJ2T = 0.0;
		for (int n = 0; n < 3; n++){
			VARJ2T = VARJ2T + 4*G*G*(0.5*Math.pow(EET.get(n), 2) + Math.pow(EET.get(n+3), 2)); 
		}
		double QTRIAL = Math.sqrt(3*VARJ2T);
		
		// 4. Check for plastic admissibility
		double PHI = QTRIAL - SIGMAY;
		IVectorRO STRES = new Array1DVector(6); 
		IVectorRO RSTAVA = EPSEN;
		//System.out.println("P: " + P);

		if (PHI / SIGMAY > TOL) {
			// Plastic step: Apply return mapping - use Newton-Raphson algorithm
			// to solve the return mapping equation (Box 7.4)
			
			//System.out.println("!!!Plastic step at Gauss Point " + "!!!");

			IFPLAS = true;
			double EPBAR = EPBARN;
			for (int n = 0; n < 10; n++){
				double DENOM = -3.0*G - gp.getMaterial().DPLFUN(sigma_y0, H, EPBAR);;
				double DDGAMA = -PHI/DENOM;
				DGAMA = DGAMA + DDGAMA;
				
				EPBAR = EPBAR + DDGAMA;
				SIGMAY = gp.getMaterial().PLFUN(sigma_y0, H, EPBAR);
				PHI = QTRIAL - 3.0*G*DGAMA - SIGMAY;
				
				double RESNOR = Math.abs(PHI / SIGMAY);
				if (RESNOR <= TOL){
					
					// update stress components
					double FACTOR = 2.0*G * (1.0 - 3.0*G*DGAMA/QTRIAL);
					STRES = EET.multiply(FACTOR).add(SOID.multiply(P));
					//System.out.println("sigma:" + MatrixFormat.format(STRES));

					// compute converged elastic (engineering) strain components
					FACTOR = FACTOR / (2.0*G);
					RSTAVA = EET.multiply(FACTOR).add(SOID.multiply(EEVD3));
					for (int k = 0; k < 3; k++){
						((Array1DVector) RSTAVA).set(k+3,EET.multiply(FACTOR).get(k+3)*2.0);
					}
					gp.getMaterial().setSigma(STRES);
					result.add(STRES);		
					
					// update accumulated plastic strain
					gp.getMaterial().setEpsPlastic(EPBAR);
					result.add(EPBAR);
					
					break;
				}
			}
			q_plot = QTRIAL - 3.0*G*DGAMA;
		}
		else {
			// Elastic step: Update stress using linear elastic law
			STRES = EET.multiply(2.0*G).add(SOID.multiply(P));
			gp.getMaterial().setSigma(STRES);
			result.add(STRES);

			// elastic engineering strain
			double EPBAR = EPBARN;
			gp.getMaterial().setEpsPlastic(EPBAR);
			result.add(EPBAR);

			q_plot = QTRIAL;
			RSTAVA = STRAT;


			//System.out.println("sigma:" + MatrixFormat.format(STRES));
		}

		LALGVA.set(0, IFPLAS);
		q_list.add(q_plot);
		eps_plot_list.add(EPSN1.get(2));
		
		result.add(RSTAVA);
		gp.getMaterial().setEpsElastic(RSTAVA);
		return result;
	}	
	
	public void tangentUpdate (double EPBAR) {
		
		IVectorRO STRES = gp.getMaterial().getSigma();
		//double DGAMA = gp.getDGAMA();
		IMatrixRO DMATX = new Array2DMatrix(6,6);
		boolean IFPLAS = LALGVA.get(0);

		if (IFPLAS == true){
			// Compute elastoplastic consistent tangent
			double R003D2 = Math.sqrt(3.0/2.0);
			// Hydrostatic pressure
			double P = 1.0/3.0 * (STRES.get(0) + STRES.get(1) + STRES.get(2));	
			// Deviatoric stress components
			IVectorRO S = STRES.subtract(SOID.multiply(P));				
			
			// Recover last elastic trial von Mises effective stress
			double SNORM = 0.0;
			for (int j = 0; j <3; j++){
				SNORM = SNORM + Math.pow(S.get(j),2) + 2.0*Math.pow(S.get(j+3), 2);
			}
			SNORM = Math.sqrt(SNORM);
			double Q = R003D2*SNORM;
			
			double QTRIAL = Q + 3.0*G*DGAMA;
			
			// Assemble elastoplastic tangent (upper triangle only)
			double AFACT = 2.0*G*(1.0 - 3.0*G*DGAMA/QTRIAL);
			double BFACT = 6*G*G*(DGAMA/QTRIAL - 1.0/(3.0*G + gp.getMaterial().DPLFUN(sigma_y0, H, EPBAR)))/(SNORM*SNORM);

			for (int m = 0; m < 6; m++){
				for (int n = m; n < 6; n++){
					((Array2DMatrix) DMATX).set(m,n,DEVPRJ.get(m, n)*AFACT 
							+ S.get(m)*S.get(n)*BFACT + K*SOID.get(m)*SOID.get(n));
				}
			}
		}
		else {
			// Compute elasticity matrix (upper triangle only)
			for (int m = 0; m < 6; m++){
				for (int n = 0; n < 6; n++){
					((Array2DMatrix) DMATX).set(m,n,DEVPRJ.get(m, n)*2.0*G 
							+ K*SOID.get(m)*SOID.get(n));
				}
			}
		}
		
		// Assemble lower triangle
		for (int n = 0; n < 5; n++){
			for (int m = n+1; m < 6; m++){
				((Array2DMatrix) DMATX).set(m,n,DMATX.get(n, m));
			}
		}
		if (IFPLAS == true){
			//System.out.println("DMATX:\n" + MatrixFormat.format(DMATX));
		}
		gp.getMaterial().setTangent(DMATX);
		
	}		

	public PlotPanel drawPlot(){
		List<Double> e = eps_plot_list;
		List<Double> s = q_list;
		DynamicDataset data = new DynamicDataset();

		for (int i = 0; i < e.size(); i++){
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