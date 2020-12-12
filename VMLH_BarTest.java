package models;

import inf.v3d.view.Viewer;

import java.util.ArrayList;
import java.util.List;

import fem_1D.Constraint;
import fem_1D.Force;
import fem_1D.GaussPoint;
import fem_1D.Node;
import plasticity.VMLH3D;
import shape_functions.LinearSF;
import shape_functions.ShapeFunctions;
import solver.NewtonRaphsonPla;
import fem_3D.Brick;
import fem_3D.BrickStructure;
import fem_3D.Visualizer3D;
import fem_abstract.Element;

public class VMLH_BarTest {
public static void main (String args[]){
		
		BrickStructure bstr = new BrickStructure();
		ShapeFunctions sf = new LinearSF("Brick");
		bstr.setSF(sf);

		for (int i = 0; i <= 10; i++){
			bstr.addNode(0.0,0.0,1.0*i);
			bstr.addNode(1.0,0.0,1.0*i);
			bstr.addNode(1.0,1.0,1.0*i);
			bstr.addNode(0.0,1.0,1.0*i);
		}
		
		Constraint c1 = new Constraint(false, false, false);
		
		for (int i = 0; i < 4; i++){
			bstr.getLON().get(i).setConstraint(c1);
			bstr.getLON().get(i+40).setConstraint(c1);
		}
		
		Force f = new Force (3.0, 0.0, 0.0);
		
		double E = 250; double nuy = 0.2;
		double yield_stress = 20;
		double hardening = 3*E/5.0;
		
		for (int i = 0; i < 10; i++){
			List<Node> ele = new ArrayList<Node>();
			for (int j = 0; j < 8; j++){
				ele.add(bstr.getLON().get(j+4*i));
			}
			Brick brick = new Brick (E, 0.3333, ele);
			brick.assignGaussPoints(2, 2, 2, sf);
			bstr.addElement(brick);
		}
		
		/*List<Node> ele1 = new ArrayList<Node>();
		List<Node> ele2 = new ArrayList<Node>();
		List<Node> ele3 = new ArrayList<Node>();
		for (int i = 0; i < 8; i++){
			ele1.add(bstr.getLON().get(i));
			ele2.add(bstr.getLON().get(i+4));
			ele3.add(bstr.getLON().get(i+8));
		}
		
		/*Brick brick1 = new Brick (E, nuy, ele1);
		brick1.assignGaussPoints(2, 2, 2, sf);
		Brick brick2 = new Brick (E, nuy, ele2);
		brick2.assignGaussPoints(2, 2, 2, sf);
		Brick brick3 = new Brick (E, nuy, ele3);
		brick3.assignGaussPoints(2, 2, 2, sf);
		
		bstr.addElement(brick1);
		bstr.addElement(brick2);
		bstr.addElement(brick3);*/
		
		bstr.assignStartIndex();
		
		for (Element b : bstr.getLOE()){
			for (GaussPoint gp : b.getLOGP()){
				VMLH3D vm = new VMLH3D(gp,yield_stress,hardening,E,nuy);
				gp.setPlasticModel(vm);
			}
		}
		
		int m = bstr.getLON().size();
		bstr.getLON().get(m/2-1).setForce(f);
		bstr.getLON().get(m/2-5).setForce(f);
		
		bstr.computeBMatrix(sf);
		
		System.out.println("************************************************************");
		System.out.println("*********************VON-MISES MATERIAL*********************");
		System.out.println("************************************************************");
		System.out.println();
				
		NewtonRaphsonPla returnMapping = new NewtonRaphsonPla(bstr);
		returnMapping.solve(5, 15.0/20.0, 1.0/20.0);
		
		/*PlotPanel plot = vm1.drawPlot();
		JFrame frame = new JFrame();	    	    
	    frame.add(plot);
	    frame.setTitle("Stress strain curve - Sigma_y = " + yield_stress + 
	    		" - Hardening = " + hardening + " - Step = " + load_step.size() );
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(800,600);
        frame.setLocation(1000,200);
        //frame.setVisible(true);*/
        
        Viewer viewer = new Viewer ();
		Visualizer3D viz = new Visualizer3D (bstr,viewer);
		viz.setConstraintSymbolScale (0.08);
		viz.setForceSymbolScale (0.5);
		viz.setForceSymbolRadius (0.04);
		viz.setDisplacementScale (1);
		viz.drawElements ();
		viz.drawConstraints();
		viz.drawForces();
		viz.drawDisplacements();
		viewer.setVisible(true);
	}
}