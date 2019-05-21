package jmetal.problems;

import java.awt.List;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.concurrent.Executor;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.encodings.solutionType.BinarySolutionType;
import jmetal.encodings.solutionType.DiscreteSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.init.Params;
import jmetal.util.JMException;

import weka.classifiers.Classifier;
import weka.classifiers.functions.MultilayerPerceptron;
import weka.classifiers.pmml.consumer.NeuralNetwork;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

public class ExternEvaluatorNN extends Problem {

	int numOfFeatures;
	
	public ExternEvaluatorNN() throws ClassNotFoundException
	{
		numberOfVariables_   = Params.getNumberOfVariables() ;
		numberOfFeatures_   = Params.getNumberOfFeatures() ;
	    numberOfObjectives_  = Params.getNumberOfObjectives() ;
	    //numberOfConstraints_ = 11;
	    problemName_         = "ExternEvaluator";
	    upperLimit_ = Params.getVariablesUpperLimit();		//new double[numberOfVariables_];
	    lowerLimit_ = Params.getVariablesLowerLimit();		//new double[numberOfVariables_];
	    
	    solutionType_ = new DiscreteSolutionType(this) ;
    	//solutionType_ = new RealSolutionType(this) ;
	    
	    numOfFeatures = Params.getNumberOfFeatures();
	}
	
	public void evaluateConstraints(Solution solution) throws JMException {
		
		return;
		//double total = 0;
	    //int number = 0;	    
	        
	    //solution.setOverallConstraintViolation(total);    
	    //solution.setNumberOfViolatedConstraint(number);  
	
	}
	
	@Override
	public void evaluate(Solution solution) throws JMException {
		List loListOfResults = new List();
        double ldRezult = 0;
        int liFeasible = 1;
        double ldDistanceToFisibility = 0;
        List loListOfFeatures = new List();
        String lsSimInDirectory;
        String lsSimOutDirectory;
        String lsSimOutFileName;
        FileOutputStream fos;
        OutputStreamWriter osw;
        BufferedWriter bw;
        FileInputStream fileInputStream;
        InputStreamReader inputStreamReader;
        BufferedReader bufferedReader;
        File simIn = new File(Params.getSimulatorSimIn());
        File simOut = new File(Params.getSimulatorSimOut());
        
        lsSimInDirectory = simIn.getParent();
        lsSimOutDirectory = simOut.getParent();         
        lsSimOutFileName = Params.getSimulatorSimOut();
        String lsSimmulatorDirectory = new File(Params.getSimulatorCommand()).getParent(); 
        File simmulatorDirectory = new File(lsSimmulatorDirectory);
        
        
        if (simIn.getParentFile().exists() == false)
            simIn.getParentFile().mkdir();

        if (simOut.getParentFile().exists() == false)
            simOut.getParentFile().mkdir();

        
      //Write parameters to simIm file
        try {
	        fos   = new FileOutputStream(Params.getSimulatorSimIn(),false);
	        osw = new OutputStreamWriter(fos);
	        bw      = new BufferedWriter(osw);
	        
        
			for (int i = 0; i < Params.getNumberOfVariables(); i++)
			{
				bw.write(solution.getDecisionVariables()[i].toString());
				bw.newLine();
			}
			/* Close the file */
			bw.close();
        }
		catch (IOException e) {
			e.printStackTrace();
		}
        
       
        //run the simulator
        try {
        	ProcessBuilder pp = new ProcessBuilder(Params.getSimulatorCommand());
        	pp.directory(simmulatorDirectory);
        	Process process = pp.start();
        	process.waitFor();
        	
		} 
        catch (Exception e) {
        	System.out.println("Command file doesn't exist.");
        
			e.printStackTrace();
		} 
        
        try
        {
        	fileInputStream = new FileInputStream(Params.getSimulatorSimOut());
        	inputStreamReader = new InputStreamReader(fileInputStream);
        	bufferedReader      = new BufferedReader(inputStreamReader); 
        	
        	liFeasible = Integer.parseInt(bufferedReader.readLine());
        	ldDistanceToFisibility =  -1 * Double.parseDouble(bufferedReader.readLine());	//-1 because OveralConstraintViolation must be negative 
        	if(ldDistanceToFisibility != 0)
        	{
        		solution.setNumberOfViolatedConstraint(1);
        		solution.setOverallConstraintViolation(ldDistanceToFisibility);
        	}
        	for(int i=0; i< Params.getNumberOfObjectives(); i++)
        		solution.setObjective(i, Double.parseDouble(bufferedReader.readLine()));
        	
        	for(int i=0; i< Params.getNumberOfFeatures(); i++)
        		solution.setFeature(i, Double.parseDouble(bufferedReader.readLine()));
	        
        	bufferedReader.close();
        }
        
        catch (Exception e) {
        	for(int i=0; i< Params.getNumberOfObjectives(); i++)
        		solution.setObjective(i, -999);
        	
        	for(int i=0; i< Params.getNumberOfFeatures(); i++)
        		solution.setFeature(i, -999);
		}
        
	}


	public void evaluateApproximateValue(Object o, Solution solution) throws JMException {
		
		MultilayerPerceptron [] tableOfNeuralNetworks = (MultilayerPerceptron []) o;
		Instance iExample;
		double evaluatedValue;
		double ldDistanceToFisibility;

		for(int i=0; i<Params.getNumberOfObjectives(); i++)
      	{      		
      		try
            {           			      			      			
      			iExample = new Instance(Params.getNumberOfVariables() + 1);
      			for(int k=0; k<Params.getNumberOfVariables();k++)      				
      				iExample.setValue(k, solution.getDecisionVariables()[k].getValue());
      					
      			//iExample.setValue(Params.getNumberOfVariables(), solution.getObjective(i));
      			//iExample.setValue(Params.getNumberOfVariables(), 99);
      			       	
      			iExample.setDataset(getInitializedInstances());
      			
      			// evaluate
      			evaluatedValue = tableOfNeuralNetworks[i].classifyInstance(iExample);
      			solution.setObjective(i, evaluatedValue);  
      			solution.setFeature(i, -1);
                    	            	
      		} 
            catch (Exception e) 
            { 
            	e.printStackTrace(); 	
            }       
      	}
		
		ldDistanceToFisibility = -1 * calulateDistanceToFeasibility(solution);
		if(ldDistanceToFisibility != 0)
    	{
    		solution.setNumberOfViolatedConstraint(1);
    		solution.setOverallConstraintViolation(ldDistanceToFisibility);
    	}
	}
	
	private double calulateDistanceToFeasibility(Solution solution) {

		double ldDistanceToFisibility =0;
		double value1;
		double value2;
		double value3;
		
		value1 = solution.getObjective(0);
		value2 = solution.getObjective(1);
		value3 = solution.getObjective(2);
		
		//first objective
		if(value1 > 1)		
			ldDistanceToFisibility = ldDistanceToFisibility + ((value1 -1)/1);		
		//else if(value1 > 11)
		//	ldDistanceToFisibility = ldDistanceToFisibility + ((value1 - 11)/1);
			
		//first objective
		if(value2 > 4)		
			ldDistanceToFisibility = ldDistanceToFisibility + ((value2 - 4)/4);		
		//else if(value2 > 15)
			//ldDistanceToFisibility = ldDistanceToFisibility + ((value2 - 15)/4);
				
		//first objective
		if(value3 > 7.5)		
			ldDistanceToFisibility = ldDistanceToFisibility + ((value3 - 7.5)/15);		
		//else if(value3 > 1130)
			//ldDistanceToFisibility = ldDistanceToFisibility + ((value3 - 1130)/15);
		
		return ldDistanceToFisibility;
	}

	private Instances getInitializedInstances() {

        FastVector attr = new FastVector(Params.getNumberOfVariables() + 1);	//+1 for the criterion
        Attribute Attribute1;
        for(int i=0; i< Params.getNumberOfVariables(); i++)
        {
        	Attribute1 = new Attribute(i+ "- input variable");
        	attr.addElement(Attribute1);
        }
        
        Attribute1 = new Attribute("- class");
    	attr.addElement(Attribute1);
    	
    	
		
		Instances iExamples = new Instances("Examples", attr, 0); 
		iExamples.setClassIndex(Params.getNumberOfVariables());
		return iExamples;
	}
}
