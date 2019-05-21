package jmetal.problems;

import java.awt.List;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
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

public class StoreSteel extends Problem {

	Hashtable<String, double[]> tableOfSolutions;	
	int accuracy = 100000;
	
	public StoreSteel() throws ClassNotFoundException
	{
		  numberOfVariables_   = 4;		
		  numberOfFeatures_   = 3;
	    numberOfObjectives_  = 3;
	    //numberOfConstraints_ = 11;
	    problemName_         = "StoreSteel";
	    upperLimit_ = new double[]{1.99, 35, 40, 65};		//Params.getVariablesUpperLimit();		//new double[numberOfVariables_];
	    lowerLimit_ = new double[]{1.5, 33, 10, 25};		//new double[numberOfVariables_];
	    
	    solutionType_ = new DiscreteSolutionType(this) ;
    	//solutionType_ = new RealSolutionType(this) ;
	    
	    Params.setDiscretizationStep(new double[]{0.01, 1, 5, 5});
	    
	    init();
	}
	
	private void init() {
	  
		BufferedReader br;
		Scanner s;
		String fileName = "ES all.txt";
		String line;
		String value;
		double doubleValue = 0;
		int intValue;		
		String variables; 
		double[] features = new double[numberOfFeatures_];
		int VariablesHash = 0;
		
		tableOfSolutions = new Hashtable<>();
		
		try {		
	    br=new BufferedReader(new FileReader(fileName));
	    while((line=br.readLine())!=null){
	    	s=new Scanner(line);
	    	variables = "";
       	for(int i=0; i< numberOfVariables_; i++)
       	{
       		value = s.next();
       		doubleValue = Double.parseDouble(value);       
       		intValue = (int) (doubleValue * accuracy);
       		variables = variables + intValue + "|";
       		VariablesHash = variables.hashCode();
       	}	  
       	
       	for(int i=0; i< numberOfFeatures_; i++)
       	{
       		value = s.next();
       		doubleValue = Double.parseDouble(value);       
       		features [i] = doubleValue;
       	}	  
       	
       	tableOfSolutions.put(variables, new double[]{features[0], features[1], features[2]});     
	    }
	    
    } catch (IOException e) {
	    
	    e.printStackTrace();
    } 
	  
  }

	public void evaluateConstraints(Solution solution) throws JMException {
		
		return;
		
	}
	
	@Override
	public void evaluate(Solution solution) throws JMException {
		
		double[] features;
		double doubleValue;
		int intValue;
		String variables = ""; 
		double distanceToFisibility;
		
		for(int i=0; i<numberOfVariables_; i++)
		{
   		doubleValue = solution.getDecisionVariables()[i].getValue();       
   		intValue = (int) (doubleValue * accuracy);
   		variables = variables + intValue + "|";			
		}
		
		features = tableOfSolutions.get(variables);
		
		for(int i=0; i<numberOfFeatures_; i++)
		{
			solution.setFeature(i, features[i]);
			solution.setObjective(i, getObjective(i, features[i]));
		}
		
		distanceToFisibility =  calculateDistanceToFeasibility(solution);
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

	private double calculateDistanceToFeasibility(Solution solution) {

		double distance = 0;
		
		if(solution.getFeature(0) < 10)
			distance = distance + 10 - solution.getFeature(0);
		else if(solution.getFeature(0) > 11)
			distance = distance + solution.getFeature(0) - 11;
		
		if(solution.getFeature(1) < 11)
			distance = distance + 11 - solution.getFeature(1);
		else if(solution.getFeature(1) > 15)
			distance = distance + (solution.getFeature(1) - 15)/4;
		
		if(solution.getFeature(2) < 1115)
			distance = distance + 1115 - solution.getFeature(2);
		else if(solution.getFeature(2) > 1130)
			distance = distance + (solution.getFeature(2) - 1130)/15;
		
	  return (-1 * distance);
  }

	private double getObjective(int i, double feature) {
		double objective;
	  
	  if(i== 0)
	  	objective = Math.abs(feature-10);
	  else if(i== 1)
	  	objective = Math.abs(feature-13);
	  else
	  	objective = Math.abs(feature-1122.5);
	  
	  return objective;
  }

}
