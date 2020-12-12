package solver;

import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVectorRO;
import iceb.jnumerics.MatrixFormat;

import java.util.ArrayList;
import java.util.List;

import fem_1D.GaussPoint;
import linalg.Array1DVector;
import fem_abstract.Element;
import fem_abstract.Structure;

public class NewtonRaphsonPla {
	private Structure str;
	
	public NewtonRaphsonPla (Structure str){
		this.str = str;
	}
	
	@SuppressWarnings("rawtypes")
	public void solve(int steps, double lambda, double delta_lambda) {
		
		IVectorRO epsConverged = null;
		IVectorRO sigConverged = null;
		IVectorRO EPSEN = null;
		double EPBARN = 0.0; 
				
		if (str.getClass().getName() == "fem_3D.BrickStructure"){
			epsConverged = new Array1DVector(6);
			sigConverged = new Array1DVector(6);
			EPSEN = new Array1DVector(6);
		}
		
		if (str.getClass().getName() == "fem_2D.QuadStructure"){
			epsConverged = new Array1DVector(4);
			sigConverged = new Array1DVector(4);
			EPSEN = new Array1DVector(4);
		}
		
		IVectorRO u_0 = new Array1DVector(str.enumerateDOFs());
		
		for (int j = 0; j < steps; j++){
			
			int count = j+1;
			System.out.println("**************************** LOAD STEP " + count + " *****************************");
			
			lambda = lambda + delta_lambda;
			IVectorRO r_ext = str.assembleExternalForces().multiply(lambda);
			 
			IVectorRO u_n = new Array1DVector(str.enumerateDOFs());
			u_n = u_0;
			
			IVectorRO u_k = u_n;
			List result = new ArrayList(); 

			double eta = 1.0;

			for (int k=0; k < 10; ++k) {
			
				IVectorRO delta_u = new Array1DVector(str.enumerateDOFs());

				if (eta > 1.0e-6){
					System.out.println();
					System.out.println("STEP " + (k+1) ); 
										
					str.setDisplacement();
					
					for (Element q : str.getLOE()){
						for (GaussPoint gp : q.getLOGP()){
							
							IVectorRO epsk = gp.getBMatrix().multiply(q.getDisplacement());
							IVectorRO DEPS = epsk.subtract(epsConverged);
//							System.out.println("epsk " + MatrixFormat.format(epsk));
							result = gp.getPlasticModel().returnMapping(epsk, DEPS, EPSEN, EPBARN, sigConverged);	
							gp.getPlasticModel().tangentUpdate((double) result.get(2));
							
						}
					}
			
					IMatrixRO K_t = str.assemblePlasticStiffnessMatrix();
					IVectorRO r_int = str.assembleInternalForces();
					
//					System.out.println("r_int:" + MatrixFormat.format(r_int));
//					System.out.println("K_t:" + MatrixFormat.format(K_t));

					LUDecomposition lud = new LUDecomposition(K_t);
					delta_u = lud.solveFor(r_ext.subtract(r_int));
					
					str.selectPlasticDisplacements(delta_u);

					u_k = u_k.add(delta_u);
//					System.out.println("delta_u " + MatrixFormat.format(delta_u) + "\n");

					eta = delta_u.normTwo()/((u_k.subtract(u_n)).normTwo());
					System.out.println("eta " + eta);
					//System.out.println();
					
				}
				
				else {
					System.out.println("RETURN MAPPING converged after " + k + " steps");
					System.out.println();
					break;
				}
				System.out.println("---------------------------------------------------------");

			}
			u_0 = u_k;
			eta = 100;
			epsConverged = (IVectorRO) result.get(0);
			sigConverged = (IVectorRO) result.get(1);
			EPBARN = (double) result.get(2);
			EPSEN = (IVectorRO) result.get(3);

			System.out.println("u_0 " + MatrixFormat.format(u_0));

			if (j == steps -1) {
				for (Element q : str.getLOE()){
					for (GaussPoint gp : q.getLOGP()){
						System.out.println("EpsElastic: " + MatrixFormat.format( gp.getMaterial().getEpsElastic()) ) ;
						System.out.println("EpsPlastic: " + ( gp.getMaterial().getEpsPlastic()) ) ;
						System.out.println("Sigma: " + MatrixFormat.format(gp.getMaterial().getSigma()) );
						System.out.println();
					}
				}	
			}
		}
	}
	
}