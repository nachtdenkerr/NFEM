package models;

import java.util.ArrayList;
import java.util.List;

import fem_1D.Constraint;
import fem_1D.Force;
import fem_1D.GaussPoint;
import fem_1D.Material;
import fem_1D.Node;
import plasticity.VMLH3D;
import shape_functions.LinearSF;
import shape_functions.ShapeFunctions;
import solver.NewtonRaphsonPla;
import fem_3D.Brick;
import fem_3D.BrickStructure;
import fem_abstract.Element;
import fem_abstract.Structure;

public class VMLH_BrickTest {
	public static void main (String args[]){
		
		Structure bstr = new BrickStructure();
		ShapeFunctions sf = new LinearSF("Brick");
		bstr.setSF(sf);
		
		Node n1 = bstr.addNode(0.0,0.0,0.0);
		Node n2 = bstr.addNode(1.0,0.0,0.0);
		Node n3 = bstr.addNode(1.0,1.0,0.0);
		Node n4 = bstr.addNode(0.0,1.0,0.0);
		Node n5 = bstr.addNode(0.0,0.0,1.0);
		Node n6 = bstr.addNode(1.0,0.0,1.0);
		Node n7 = bstr.addNode(1.0,1.0,1.0);
		Node n8 = bstr.addNode(0.0,1.0,1.0);
		List<Node> l = new ArrayList<Node>();
		l.add(n1); l.add(n2); l.add(n3); l.add(n4); l.add(n5); l.add(n6); l.add(n7); l.add(n8);
		
		Constraint c1 = new Constraint(false, false, false);
		Constraint c2 = new Constraint(true,  false, false);
		Constraint c3 = new Constraint(true,  true,  false);
		Constraint c4 = new Constraint(false, true,  false);
		Constraint c5 = new Constraint(false, false, true );
		Constraint c6 = new Constraint(true,  false, true );
		Constraint c7 = new Constraint(true,  true,  true );
		Constraint c8 = new Constraint(false, true,  true );
		
		n1.setConstraint(c1); n2.setConstraint(c2); n3.setConstraint(c3); n4.setConstraint(c4);
		n5.setConstraint(c5); n6.setConstraint(c6); n7.setConstraint(c7); n8.setConstraint(c8);
		
		Force f = new Force (0.0, 0.0, 8.0);
		n5.setForce(f);
		n6.setForce(f);
		n7.setForce(f);
		n8.setForce(f);
				
		double E = 250, nuy = 0.3333;
		double yield_stress = 20.0;
		double hardening = 1*E/5.0;
		
		Element brick = new Brick (E, 0.3333, l);
		((Brick) brick).assignGaussPoints(2, 2, 2, sf);
		
		bstr.addElement(brick);
		bstr.assignStartIndex();
		
		for (Element b : bstr.getLOE()){
			for (GaussPoint gp : b.getLOGP()){
				Material mat = new Material(1);
				gp.setMaterial(mat);
				VMLH3D vm = new VMLH3D(gp,yield_stress,hardening,E,nuy);
				gp.setPlasticModel(vm);
			}
		}
		
		brick.computeBMatrix(sf);
		
		NewtonRaphsonPla returnMapping = new NewtonRaphsonPla(bstr);
		returnMapping.solve(8, 1.0/8.0, 1.0/8.0);
		
		/*PlotPanel plot = vm.drawPlot();
		JFrame frame = new JFrame();	    	    
	    frame.add(plot);
	    frame.setTitle("VMLH_Brick_Test - "+"Stress strain curve - Sigma_y = " + yield_stress + 
	    		" - Hardening = " + hardening + " - Step = " + load_step.size() );
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(800,600);
        frame.setLocation(1000,200);
        //frame.setVisible(true);
        
        Viewer viewer = new Viewer ();
		Visualizer3D viz = new Visualizer3D (bstr,viewer);
		viz.setConstraintSymbolScale (0.08);
		viz.setForceSymbolScale (0.05);
		viz.setForceSymbolRadius (0.008);
		viz.setDisplacementScale (1);
		viz.drawElements ();
		viz.drawConstraints();
		viz.drawForces();
		viz.drawDisplacements();*/
		//viewer.setVisible(true);
	}
}