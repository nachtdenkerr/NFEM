package models;

import fem_1D.Constraint;
import fem_1D.Force;
import fem_1D.GaussPoint;
import fem_1D.Material;
import fem_1D.Node;
import fem_2D.Quad;
import fem_2D.QuadStructure;
import fem_abstract.Element;
import plasticity.VMLH2D;
import solver.NewtonRaphsonPla;

public class VMLH_QuadTest {
	public static void main (String args[]){
		Node n1 = new Node (0.0, 0.0, 0.0);
		Node n2 = new Node (1.0, 0.0, 0.0);
		Node n3 = new Node (1.0, 1.0, 0.0);
		Node n4 = new Node (0.0, 1.0, 0.0);
		
		Constraint c1 = new Constraint(false, false, false);
		Constraint c2 = new Constraint(true , false, false);
		Constraint c3 = new Constraint(true , true , false);
		Constraint c4 = new Constraint(false, true , false);
		
		n1.setConstraint(c1);
		n2.setConstraint(c2);
		n3.setConstraint(c3);
		n4.setConstraint(c4);
		
		double E = 210.0e6, nuy = 0.3;
		Quad quad = new Quad(E, 0.3333, "plane strain", 1.0, n1, n2, n3, n4);
		
		Force f = new Force (0.0, 135000.0, 0.0);
		n3.setForce(f);
		n4.setForce(f);
			
		QuadStructure qstr = new QuadStructure();
		qstr.addElement(quad);
		qstr.assignStartIndex();
		//qstr.computeB_Plastic_Matrix();
		
		double yield_stress = 235.0e3;
		double hardening = 21000.0;

		qstr.assignPlasticGaussPoints(4);

		for (Element q : qstr.getLOE()) {
			for (GaussPoint gp : q.getLOGP()){
				Material mat = new Material(1);
				gp.setMaterial(mat);
				VMLH2D vm = new VMLH2D(gp, yield_stress, hardening, E, nuy);
				gp.setPlasticModel(vm);
			}
		}
		
		
		
		System.out.println("************************************************************");
		System.out.println("*********************VON-MISES MATERIAL*********************");
		System.out.println("************************************************************");
		System.out.println();
		
		NewtonRaphsonPla returnMapping = new NewtonRaphsonPla(qstr);
		returnMapping.solve(10, 0, 1);
		
		
		
		/*PlotPanel plot = qstr.getLOE().get(0).getLOGP().get(0).getPlasticModel().drawPlot();
		JFrame frame = new JFrame();	    	    
	    frame.add(plot);
	    frame.setTitle("Stress strain curve - Sigma_y = " + yield_stress + 
	    		" - Hardening = " + hardening + " - Step = " );
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(800,600);
        frame.setLocation(1000,200);*/
        //frame.setVisible(true);
        
        /*Viewer viewer = new Viewer ();
		Visualizer2D viz = new Visualizer2D (qstr,viewer);
		viz.setConstraintSymbolScale (0.08);
		viz.setForceSymbolScale (0.05);
		viz.setForceSymbolRadius (0.008);
		viz.setDisplacementScale (100);
		viz.drawElements ();
		viz.drawConstraints();
		viz.drawForces();
		viz.drawDisplacements();
		//viewer.setVisible(true);*/
	}
}