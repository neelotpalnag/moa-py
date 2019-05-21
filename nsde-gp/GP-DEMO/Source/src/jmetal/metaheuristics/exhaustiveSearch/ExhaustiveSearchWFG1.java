/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jmetal.metaheuristics.exhaustiveSearch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Properties;

import jmetal.core.Variable;
import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.NonDominatedSolutionList;
import jmetal.core.Problem;
import jmetal.core.Algorithm;
/**
 *
 * @author Rok
 */
public class ExhaustiveSearchWFG1 extends Algorithm {
   
	String logFileName;
  /**
  * Constructor
  * @param problem Problem to solve
  */
  public ExhaustiveSearchWFG1(Problem problem){
    super (problem) ;
  } // ExhaustiveSearch
  
  /**
   * Constructor 
   * @param problem Problem to solve
   * @param properties Properties of algorithm
   */
  public ExhaustiveSearchWFG1(Problem problem, Properties properties){
    super (problem, properties);  
  }
  
   /**
  * Runs the ExhaustiveSearch algorithm.
  * @return a <code>SolutionSet</code> that is a set of solutions
  * as a result of the algorithm execution
   * @throws JMException
  */
  
    @Override
   public SolutionSet execute() throws JMException, ClassNotFoundException {        
        
        NonDominatedSolutionList ndl = new NonDominatedSolutionList();

        Variable[] variables = new Variable[5];   
        Solution solution = new Solution(problem_);
        
        variables = solution.getDecisionVariables();
        variables[2].setValue(2.1);
        variables[3].setValue(2.8);
        variables[4].setValue(3.5);
        variables[5].setValue(4.2);
        
        logFileName = getParameter("Output.logFileName");
        
        deleteDejan();  
        
        for(double prvi=0; prvi<=0.00001; prvi=prvi+0.000000001)
        {
        	for(double drugi = 2; drugi <= 2.2; drugi = drugi + 0.1)
        	{
        		variables[0].setValue(prvi);
        		variables[1].setValue(drugi);
        		solution.setDecisionVariables(variables);
        		problem_.evaluate(solution);        		
        		ndl.add(new Solution(solution));
        	}
        }
        
        for(double prvi=0.1; prvi<=2; prvi=prvi+0.2)
        {
        	for(double drugi = 0; drugi <= 4; drugi = drugi + 0.2)
        	{
        		variables[0].setValue(prvi);
        		variables[1].setValue(drugi);
        		solution.setDecisionVariables(variables);
        		problem_.evaluate(solution);        		
        		ndl.add(new Solution(solution));
        	}
        }
        
        ndl.printVariablesAndObjectivesToFile(logFileName + "dejan.txt");
        
        
        
        
        
        
        
        return ndl;
   }//execute
    
    private void deleteDejan() {
    	try {
      	BufferedWriter bw;    
  	    bw = new BufferedWriter(new FileWriter(logFileName + "dejan.txt"));    
      bw.write("");
      bw.flush();
      bw.close();
      } catch (IOException e) {
  	    e.printStackTrace();
      }
    }
}//RandomSearch
