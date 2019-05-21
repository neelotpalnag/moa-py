package jmetal.problems;

import java.awt.List;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.lang.ProcessBuilder.Redirect;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Scanner;
import java.util.concurrent.Executor;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.encodings.solutionType.BinarySolutionType;
import jmetal.encodings.solutionType.DiscreteSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.init.Params;
import jmetal.util.JMException;

import java.security.*;

public class EKGSimulator extends Problem {

	Hashtable<String, double[]> tableOfSolutions;	
	int accuracy = 100000;
	
	public EKGSimulator() throws ClassNotFoundException
	{
		  numberOfVariables_   = 8;		
		  numberOfFeatures_   = 0;
	    numberOfObjectives_  = 2;
	    //numberOfConstraints_ = 11;
	    problemName_         = "EKGSimulator";
	    lowerLimit_ = new double[]{0.0003, 0.01,  0.01, 330, 0.0003, 0.01,  0.01, 330};		//Params.getVariablesUpperLimit();		//new double[numberOfVariables_];
	    upperLimit_ = new double[]{0.0010, 0.10,  0.10, 400, 0.0010, 0.10,  0.10, 400};		//new double[numberOfVariables_];
	  		    
	    //solutionType_ = new DiscreteSolutionType(this) ;
    	solutionType_ = new RealSolutionType(this) ;
	    
	    //Params.setDiscretizationStep(new double[]{0.01, 1, 5, 5});
	    
	    init();
	}
	
	private void init() {
	  		
  }

	public void evaluateConstraints(Solution solution) throws JMException {
		
		return;
		
	}
	
	@Override
	public void evaluate(Solution solution) throws JMException {
		
		String simulatorCommand = System.getProperty("user.dir") + "\\EKGSimulator\\singleRun.bat";
		File simmulatorDirectory = new File(System.getProperty("user.dir") + "\\EKGSimulator\\");	//"D:\\IJS\\EKG_simulator_MDepolli\\");
		File outputFile = new File(System.getProperty("user.dir") + "\\EKGSimulator\\izhod.txt");
		String objectives = "";
		String violation = "";
		
		try{
		  // Create file 
		  FileWriter fstream = new FileWriter(simulatorCommand);
		  BufferedWriter out = new BufferedWriter(fstream);
		  out.write("ekgsim.exe -sim ");
		  out.write(String.valueOf(solution.getDecisionVariables()[0].getValue()));
		  for(int i=1; i< numberOfVariables_; i++)
		  	out.write("," + String.valueOf(solution.getDecisionVariables()[i].getValue()));
		  //Close the output stream
		  out.close();
		  fstream.close();
		}
		catch (Exception e){//Catch exception if any
		  System.err.println("Error: " + e.getMessage());
		}
		
		try
		{
			ProcessBuilder pp = new ProcessBuilder(simulatorCommand);
    	pp.directory(simmulatorDirectory);
    	pp.redirectOutput(outputFile);
    	Process process = pp.start();
    	process.waitFor();			
		} 
		catch (IOException | InterruptedException e) {	   
	    e.printStackTrace();
    }
		finally	{		}
		
		try{
			  FileInputStream fstream = new FileInputStream(outputFile);
			  DataInputStream in = new DataInputStream(fstream);
			  BufferedReader br = new BufferedReader(new InputStreamReader(in));
			  String strLine;
	
			  while ((strLine = br.readLine()) != null)   {
			  	if(strLine.contains("criteria ="))
			  	{
			  		objectives = strLine.substring(strLine.indexOf('<')+1, strLine.indexOf('>'));
			  		violation = strLine.substring(strLine.indexOf("violation") + 12);
			  		break;
			  	}
	
			  }
			  //Close the input stream
			  in.close();
			  fstream.close();
		}			
		catch (Exception e){//Catch exception if any
			  System.err.println("Error: " + e.getMessage());
		}
		
		String[] tableOfObjectives = objectives.split(",");
		
		for(int i=0; i<numberOfObjectives_; i++)
			solution.setObjective(i, Double.parseDouble(tableOfObjectives[i]));
		
		
		double distanceToFisibility =  Double.parseDouble(violation);
		if(distanceToFisibility != 0)
  	{
  		solution.setNumberOfViolatedConstraint(1);
  		solution.setOverallConstraintViolation(distanceToFisibility);
  	}  	
  	else
  	{
  		solution.setNumberOfViolatedConstraint(0);
  		solution.setOverallConstraintViolation(0);
  	}
	}

}
