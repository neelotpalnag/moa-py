package jmetal.init;

import jmetal.core.Model;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Variable;
import jmetal.models.SPGPHT;
import jmetal.models.SPGPm;
import jmetal.models.SPGPm2;
import jmetal.problems.WFG.WFG1;
import jmetal.util.JMException;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class WFG1TestCalculateValues {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		
		int numberOfExamples = 1001;
		int numberOfVariables = 6;
		int numberOfObjectives = 2;
		int activeSet= 400;
		SolutionSet learnExamples = new SolutionSet(numberOfExamples);
		SolutionSet testExamples = new SolutionSet(numberOfExamples);
		Variable[] variables = new Variable[6];
		Problem problem;
		Solution solution;
		
		problem = new WFG1("Real");
		Params.setSeed(0);
		
    solution = new Solution(problem);
    variables = solution.getDecisionVariables();    
    
    int maxOptIter = -200;
    
    
    Model globalModel = new SPGPHT(numberOfVariables , numberOfObjectives, activeSet);
    //((SPGPm2)globalModel).setHyperparameters(null);
    ((SPGPHT)globalModel).setMaxOptIter(maxOptIter);
    ((SPGPHT)globalModel).setWindowSize(numberOfExamples);
    
    

    //read learn examples
    BufferedReader readbuffer = new BufferedReader(new FileReader("D:\\Doktorat\\GP\\input_test.txt"));
    String strRead;

    numberOfExamples = 0;
    
    while ((strRead=readbuffer.readLine())!=null){
    	String splitarray[] = strRead.split("\t");
    	variables[0].setValue(Double.parseDouble(splitarray[0]));
    	variables[1].setValue(Double.parseDouble(splitarray[1]));
    	variables[2].setValue(Double.parseDouble(splitarray[2]));
    	variables[3].setValue(Double.parseDouble(splitarray[3]));
    	variables[4].setValue(Double.parseDouble(splitarray[4]));
    	variables[5].setValue(Double.parseDouble(splitarray[5])); 
    	solution.setDecisionVariables(variables);
    	problem.evaluate(solution);        		
    	testExamples.add(new Solution(solution));
    	numberOfExamples++;
    	//System.out.println(Double.parseDouble(splitarray[0]) + " " + Double.parseDouble(splitarray[1]));
    	}

    readbuffer.close();	
        
  
  FileWriter fos = new FileWriter("D:\\Doktorat\\GP\\target_test.csv");
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
