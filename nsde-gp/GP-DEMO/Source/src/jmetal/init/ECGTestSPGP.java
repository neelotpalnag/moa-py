package jmetal.init;

import jmetal.core.Model;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Variable;
import jmetal.models.SPGPm;
import jmetal.problems.EKGSimulator;
import jmetal.problems.WFG.WFG1;
import jmetal.util.JMException;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class ECGTestSPGP {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		
		int numberOfExamples = 1000;
		int numberOfVariables = 8;
		int numberOfObjectives = 2;
		int activeSet= 400;
		SolutionSet learnExamples = new SolutionSet(numberOfExamples);
		SolutionSet testExamples = new SolutionSet(numberOfExamples);
		Variable[] variables = new Variable[numberOfVariables];
		Problem problem;
		Solution solution;
		
		problem = new EKGSimulator();
		Params.setSeed(0);
		
    solution = new Solution(problem);
    variables = solution.getDecisionVariables();    
    
    int maxOptIter = -200;
    
    
    Model globalModel = new SPGPm(numberOfVariables , numberOfObjectives, activeSet);
    ((SPGPm)globalModel).setHyperparameters(null);
    ((SPGPm)globalModel).setMaxOptIter(maxOptIter);
    ((SPGPm)globalModel).setWindowSize(numberOfExamples);
    
    

    //read learn examples
    BufferedReader readbuffer = new BufferedReader(new FileReader("D:\\Doktorat\\GP\\ECG_input_train.txt"));
    BufferedReader readbufferTarget = new BufferedReader(new FileReader("D:\\Doktorat\\GP\\ECG_target_train.txt"));
    String strRead;
    String strReadTarget;

    while ((strRead=readbuffer.readLine())!=null){
    	strReadTarget = readbufferTarget.readLine();
    	
    	String splitarray[] = strRead.split("\t");
    	String splitarrayTarget[] = strReadTarget.split("\t");
    	
    	for(int i=0; i< numberOfVariables; i++)
    	{
	    	variables[i].setValue(Double.parseDouble(splitarray[i]));	    	
	    }
    	solution.setDecisionVariables(variables);
    	
    	for(int i=0; i< numberOfObjectives; i++)
    	{
    		solution.setObjective(i, Double.parseDouble(splitarrayTarget[i]));
    	}
    	
    	//problem.evaluate(solution);        		
    	learnExamples.add(new Solution(solution));
    	//System.out.println(Double.parseDouble(splitarray[0]) + " " + Double.parseDouble(splitarray[1]));
    	}

    readbuffer.close();	
    readbufferTarget.close();	
           
    
  //train model
    globalModel.update(learnExamples);
  
  //read test examples
    readbuffer = new BufferedReader(new FileReader("D:\\Doktorat\\GP\\ECG_input_test.txt"));


    while ((strRead=readbuffer.readLine())!=null){
    	String splitarray[] = strRead.split("\t");
    	for(int i=0; i< numberOfVariables; i++)
    	{
	    	variables[i].setValue(Double.parseDouble(splitarray[i]));	    	
	    }
    	solution.setDecisionVariables(variables);
    	//problem.evaluate(solution);        		
    	testExamples.add(new Solution(solution));
    	//System.out.println(Double.parseDouble(splitarray[0]) + " " + Double.parseDouble(splitarray[1]));
    	}

    readbuffer.close();	

    
    //approximate;
	  for(int i=0; i<numberOfExamples; i++)
	  {
	  	globalModel.evaluate(testExamples.get(i));
	  }
  
  
  FileWriter fos = new FileWriter("D:\\Doktorat\\GP\\ECG_target_test.csv");
  PrintWriter dos = new PrintWriter(fos);


  for (int i=0;i<numberOfExamples;i++)
  {
  	for(int j=0; j<numberOfVariables; j++)
  		dos.print(testExamples.get(i).getDecisionVariables()[j].getValue()+"\t");
  	
  	dos.print("|\t");
  	
  	for(int j=0; j<numberOfObjectives; j++)
  		dos.print(testExamples.get(i).getObjective(j)+"\t");
  	
  	dos.print("|\t");
  	
  	for(int j=0; j<numberOfObjectives; j++)
  		dos.print(testExamples.get(i).getstandardDeviance(j)+"\t");
  	
  	dos.println();
  }
  dos.close();
  fos.close();
  
	}
}
