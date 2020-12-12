package plasticity;

import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVectorRO;
import iceb.jnumerics.MatrixFormat;

import java.util.ArrayList;
import java.util.List;

import fem_1D.GaussPoint;
import linalg.Array1DVector;
import linalg.Array2DMatrix;

public class DruckerPrager_2D extends PlasticModel{
	private GaussPoint gp;
	private double c0, theta, K, G, H, ETA, XI, ETABAR;
	private IVectorRO SOID;
	private IMatrixRO FOID, DEVPRJ;
		
	private List<Double> eps_plot_list = new ArrayList<Double>();
	private List<Double> q_list = new ArrayList<Double>();
	private List<Boolean> LALGVA = new ArrayList<Boolean>();
	
	private int MAXRT = 50;
	private double DGAM;
	private static double TOL = 1e-8;
	
	public DruckerPrager_2D (GaussPoint gp, double cohesion , double H, double E, double nuy, double angle){
		this.gp = gp; this.theta = angle*Math.PI/180.0; this.c0 = cohesion; this.H = H;

		// initialize first history variables
		IVectorRO v0_4 = new Array1DVector(0.0, 0.0, 0.0, 0.0);
		double c0 = 0.0;

		LALGVA.add(false);
		LALGVA.add(false);
		gp.getMaterial().setSigma(v0_4);
		
		eps_plot_list.add(c0);
		q_list.add(c0);
		
		K = E/(3.0*(1.0 - 2.0*nuy));
		G = E/(2.0*(1.0 + nuy));
		
		SOID = new Array1DVector(1.0, 1.0, 0.0, 1.0);

		FOID = new Array2DMatrix(4,4);
		((Array2DMatrix) FOID).set(0, 0, 1.0);
		((Array2DMatrix) FOID).set(1, 1, 1.0);
		((Array2DMatrix) FOID).set(2, 2, 0.5);
		((Array2DMatrix) FOID).set(3, 3, 1.0);
		
		DEVPRJ = new Array2DMatrix(4, 4);
		for (int a = 0; a < 4; a++){
			for (int b = 0; b < 4; b++){
				((Array2DMatrix) DEVPRJ).set(a,b,FOID.get(a, b) - (1.0/3.0)*SOID.get(a)*SOID.get(b));
			}
		}
		
		ETA = 3.0 * Math.tan(theta) / (Math.sqrt(9.0 + 12.0 * Math.pow(Math.tan(theta), 2)));
		XI = 3.0 / (Math.sqrt(9.0 + 12.0*Math.pow(Math.tan(theta), 2)));
		ETABAR = 3.0 * Math.tan(theta) / Math.sqrt(9 + 12*Math.tan(theta) * Math.tan(theta));
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	public List<IVectorRO> returnMapping(IVectorRO EPSN1, IVectorRO DEPS, IVectorRO EPSEN, double EPBARN, IVectorRO STRESN) {
		
		//1. Get history variables
		boolean IFPLAS = LALGVA.get(0);
		boolean APEX = false, SUFAIL = false;
		LALGVA.add(APEX);
		double DGAMA = 0.0;
		double EPBAR = EPBARN;
		IVectorRO RSTAVA = EPSEN;
		
		List result = new ArrayList();
		result.add(EPSN1);
		
		//Compute elastic trial state
		//---------------------------
		//Elastic trial volumetric strain and pressure stress
		IVectorRO STRAT = EPSEN.add(DEPS);
		gp.setSTRAT(STRAT);
		
		double EETV = STRAT.get(0) + STRAT.get(1) + STRAT.get(3);
		double PT = K*EETV;
		
		//Elastic trial deviatoric stress
		double EEVD3 = EETV/3.0;
		IVectorRO STRIAL = (STRAT.subtract(SOID.multiply(EEVD3)) ).multiply(2.0*G);
		((Array1DVector) STRIAL).set(2,STRAT.get(2) * 0.5 * 2.0 * G);
		System.out.println("STRAT: "+ MatrixFormat.format(STRAT));
		System.out.println("EPBAR: " +(EPBAR));
		System.out.println("ETA: " + ETA + "    XI: " + XI);
		System.out.println("strainElasticNewTrialVol: " + (EETV));

		//Compute elastic trial stress J2 invariant and cohesion
		double VARJ2T = Math.pow(STRIAL.get(2), 2) + 
				0.5*( Math.pow(STRIAL.get(0), 2) + Math.pow(STRIAL.get(1), 2) + Math.pow(STRIAL.get(3), 2) );
		double COHE = c0 + H * EPBARN;
		
		//2.Check for plastic consistency
		double SQRJ2T = Math.sqrt(VARJ2T);
		double PHI = SQRJ2T + ETA*PT - XI*COHE;
		double RES = PHI;
		double FACTOR = -100.0;
		System.out.println("J2: " + VARJ2T);
		System.out.println("Yield function: " + RES);

		if (COHE != 0.0) RES = RES/Math.abs(COHE);
		if (RES >= TOL){
			//Plastic step: Use return mapping
			//================================
			IFPLAS = true;
			APEX = false;
			
			//Apply return mapping to smooth portion of cone - Box 8.9
			//--------------------------------------------------------
			for (int j = 0; j < MAXRT; j++){
				//Compute residual derivative
				double DENOM = - G - K*ETABAR*ETA - XI*XI*gp.getMaterial().DPLFUN(c0, H, EPBAR);
				//Compute Newton-Raphson increment and update variable DGAMA
				double DDGAMA = - PHI/DENOM;
				DGAMA = DGAMA + DDGAMA;
				//Compute new residual
				EPBAR = EPBARN + XI * DGAMA;
				COHE = gp.getMaterial().PLFUN(c0, H, EPBAR);
				double SQRJ2 = SQRJ2T - G * DGAMA;
				double P = PT - K * ETABAR * DGAMA;
				PHI = SQRJ2 + ETA * P - XI * COHE;
				//Check convergence
				double RESNOR = Math.abs(PHI);
				if (COHE != 0.0) RESNOR = RESNOR / Math.abs(COHE);
//				System.out.println(j);
				
				if (RESNOR <= TOL){
					//Check validity of return to smooth portion
					if (SQRJ2 >= 0.0){
						//results are valid, update stress components and other variables
						System.out.println("!!!Smooth portion!!!");

						if (SQRJ2T == 0.0) FACTOR = 0.0;
						else {
							FACTOR = 1.0 - G * DGAMA/SQRJ2T;
							System.out.println("FACTOR: " + FACTOR);
							System.out.println("DGAMA: " + DGAMA);
							System.out.println("SQRJ2T: " + SQRJ2T);
						}
					}
					else {
						//smooth wall return not valid - go to apex return procedure
						//Apply return mapping to APEX - Box 8.10
						//--------------------------------------- 
						//perform checks and set some variables
						System.out.println("!!!Apex portion!!!");

						APEX = true;
						//if (ETA == 0.0) 
						double ALPHA = XI/ETABAR;
						double BETA = XI/ETA;
						//Set initial guess for unknown DEPV and start iterations
						double DEPV = 0.0;
						EPBAR = EPBARN;
						COHE = c0 + H * EPBAR;
						RES = BETA * COHE - PT;
						for (int k = 0; k < MAXRT; k++){
							DENOM = ALPHA * BETA * gp.getMaterial().DPLFUN(c0, H, EPBAR) + K;
							//Compute Newton-Raphson increment and update variable DEPV
							double DDEPV = -RES/DENOM;
							DEPV = DEPV + DDEPV;
							//Compute new residual
							EPBAR = EPBARN + ALPHA * DEPV;
							COHE = gp.getMaterial().PLFUN(c0, H, EPBAR);
							P = PT - K * DEPV;
							RES = BETA * COHE - P;
							//Check convergence
							RESNOR = Math.abs(RES);
							if (COHE != 0.0) RESNOR = RESNOR/Math.abs(COHE);
							if (RESNOR <= TOL) {
								DGAMA = DEPV/ETABAR;
								FACTOR = 0.0;
								break;
							}
							SUFAIL = true;
						}
					}
					//failure of stress update procedure
					SUFAIL = true;

					//Store converged stress components and other state variables
					//-----------------------------------------------------------
					IVectorRO STRES = STRIAL.multiply(FACTOR).add(SOID.multiply(P));
//					System.out.println("pNewTrial: "+ P);
//					System.out.println("sNewTrial: "+ MatrixFormat.format(STRIAL.multiply(1.0)));
//					System.out.println("FACTOR: "+ FACTOR);
					System.out.println("STRES: "+ MatrixFormat.format(STRES));

					result.add(STRES);
					gp.getMaterial().setSigma(STRES);

					gp.getMaterial().setEpsPlastic(EPBAR);

					//update EPBAR
					result.add(EPBAR);
					//compute converged elastic (engineering) strain components
					FACTOR = FACTOR/(2.0*G);
					EEVD3 = P/(K*3.0);
					RSTAVA = STRIAL.multiply(FACTOR).add(SOID.multiply(EEVD3));
					((Array1DVector) RSTAVA).set(2, FACTOR*STRIAL.get(2)*2.0);

					break;
				}
			}
		}
		else {
			//Elastic step: update stress using linear elastic law
			//====================================================
			IVectorRO STRES = STRIAL.add(SOID.multiply(PT));
			System.out.println("STRES: "+ MatrixFormat.format(STRES));
			result.add(STRES);
			gp.getMaterial().setSigma(STRES);
			
			EPBAR = EPBARN;
			gp.getMaterial().setEpsPlastic(EPBAR);
			result.add(EPBAR);
			//elastic engineering strain
			RSTAVA = STRAT;
		}
//		System.out.println();
		result.add(RSTAVA);
		gp.getMaterial().setEpsElastic(RSTAVA);

		LALGVA.set(0, IFPLAS);
		LALGVA.set(1, SUFAIL);
		LALGVA.set(2, APEX);
		DGAM = DGAMA;
		return result;
	}
	
	@Override
	public void tangentUpdate(double EPBAR){
		double DGAMA = DGAM;
		IVectorRO STRAT = gp.getSTRAT();

		double R00T2 = Math.sqrt(2.0);
		boolean APEX = LALGVA.get(2);
		boolean EPFLAG = LALGVA.get(0);
		IMatrixRO DMATX = new Array2DMatrix(4, 4);
		if (EPFLAG == true){
			//Compute elastoplastic consistent tangent
			//========================================
			//Hardening slope
			double HSLOPE = gp.getMaterial().DPLFUN(c0, H, EPBAR);
			if (APEX == true){
				//Elastoplastic tangent consistent with apex return
				//-------------------------------------------------
				double ALPHA = XI/ETABAR;
				double BETA = XI/ETA;
				double AFACT = K*(1.0 - K/(K + ALPHA*BETA*HSLOPE));
				for (int i = 0; i < 4; i++){
					for (int j = 0; j< 4; j++){
						((Array2DMatrix) DMATX).set(i,j,AFACT*SOID.get(i)*SOID.get(j));
					}
				}

			}
			else {
				//Elastoplastic tangent consistent with smooth cone wall return
				//-------------------------------------------------------------
				//Elastic trial deviatoric (physical) strain
				double EEVD3 = (STRAT.get(0) + STRAT.get(1) + STRAT.get(3))/3.0;
				IVectorRO EETD = STRAT.subtract(SOID.multiply(EEVD3));
				((Array1DVector) EETD).set(2,STRAT.get(2)*0.5);
				double ETDNOR = Math.sqrt(Math.pow(EETD.get(0), 2) + Math.pow(EETD.get(1), 2) 
						+ Math.pow(EETD.get(3), 2) + 2.0 * Math.pow(EETD.get(2), 2));
				double EDNINV;
				IVectorRO UNIDEV = new Array1DVector(4);
				
				//Unit deviatoric flow vector
				if (ETDNOR != 0.0){
					EDNINV = 1.0/ETDNOR;
				}
				else EDNINV = 0.0;
				UNIDEV = EETD.multiply(EDNINV);
				
				//Assemble tangent
				double AUX = 1.0 / (G + K*ETA*ETABAR + XI*XI*HSLOPE);
				double AFACT = 2.0 * G * (1.0 - DGAMA/(R00T2 * ETDNOR));
				double AFACD3 = AFACT/3.0;
				double BFACT = 2.0 * G * (DGAMA/(R00T2*ETDNOR) - G*AUX);
				double CFACT = -R00T2 * G * K * AUX;
				double DFACT = K * (1.0 - K * ETA * ETABAR * AUX);
				
						//DMATX = (FOID.multiply(AFACT)).add(UNIDEV.dyadicProduct(UNIDEV).multiply(BFACT))
						//		.add((UNIDEV.dyadicProduct(SOID).multiply(ETA).add(SOID.dyadicProduct(UNIDEV).multiply(ETABAR))).multiply(CFACT))
						//				.add((SOID.dyadicProduct(SOID)).multiply(DFACT - AFACD3));
				for (int i = 0; i < 4; i++){
					for (int j = 0; j < 4; j++){
						((Array2DMatrix) DMATX).set( i, j, AFACT*FOID.get(i, j) + BFACT*UNIDEV.get(i)*UNIDEV.get(j)
								+ CFACT*( ETA*UNIDEV.get(i)*SOID.get(j) + ETABAR*SOID.get(i)*UNIDEV.get(j) )
								+ (DFACT - AFACD3)*SOID.get(i)*SOID.get(j) );
					}
				}
			}
		}
		else {
			//Compute elasticity matrix
			//=========================
			double FACTOR = K - 2.0*G/3.0;
			for (int i = 0; i < 4; i++){
				for (int j = 0; j <4; j++){
					((Array2DMatrix) DMATX).set(i, j, 2.0*G*FOID.get(i, j) + FACTOR*SOID.get(i)*SOID.get(j));
				}
			}
		}
		for (int j = 0; j < 3; j++){
			for (int i = j+1; i < 4; i++){
				((Array2DMatrix) DMATX).set(i,j,DMATX.get(j, i));
			}
		}
		
		gp.getMaterial().setTangent(DMATX);;
		System.out.println("DMATX: \n"+ MatrixFormat.format(DMATX));
		System.out.println();
	}
}