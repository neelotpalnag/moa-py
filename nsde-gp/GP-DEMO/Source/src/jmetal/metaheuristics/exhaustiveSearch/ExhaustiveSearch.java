/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jmetal.metaheuristics.exhaustiveSearch;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Properties;
import jmetal.core.*;
import jmetal.encodings.variable.Real;
import jmetal.util.JMException;
import jmetal.util.NonDominatedSolutionList;
/**
 *
 * @author Rok
 */
public class ExhaustiveSearch extends Algorithm {
    
  /**
  * Constructor
  * @param problem Problem to solve
  */
  public ExhaustiveSearch(Problem problem){
    super (problem) ;
  } // ExhaustiveSearch
  
  /**
   * Constructor 
   * @param problem Problem to solve
   * @param properties Properties of algorithm
   */
  public ExhaustiveSearch(Problem problem, Properties properties){
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
        
        Solution newSolution;
        int vectorLength;
        
        ArrayList var = new ArrayList();                                //ArrayList contains vectors of each decision variable (e.g. Kursawe problem - ArrayList contains 2 vectors)
        ArrayList<Integer> vectorLengths=new ArrayList<Integer>();      //Contains vector length of each decision variable which depends of resolution (step)
        
        NonDominatedSolutionList ndl = new NonDominatedSolutionList();
        SolutionSet allSolutions = new SolutionSet();
        
        /*
         * Tukaj bo potrebno spremenit, če bomo vnašali različne korake za različne kriterije
         */
        double step=Double.parseDouble(getParameter("Problem.step"));  //Korak
        String logFileName = getParameter("Output.logFileName");
        
        int numberOfVariables=problem_.getNumberOfVariables();          //Number of decision variables
        int i=0, numberOfVectors=1;
        int [] index = new int[numberOfVariables];
        boolean [] pogoj = new boolean[numberOfVariables];
          
        while(i!=numberOfVariables){
                        
            vectorLength=(int)Math.floor(((problem_.getUpperLimit(i)-problem_.getLowerLimit(i))/step)+1);   //Number of elements inside of interval bounds according to resolution (step)
            vectorLengths.add(vectorLength);         
            var.add(new ArrayList<Double>());
            double temp=problem_.getLowerLimit(i);  //temp is containing the lower bound for each decision variable
            
            for(int j=0; j<vectorLengths.get(i); j++){  //Initialization of all elements for each vector
                ((ArrayList)var.get(i)).add(temp);      
                temp=temp+step;
            }       
            i++;
        }   
        
        Iterator iter = vectorLengths.iterator();
        while(iter.hasNext()){
            numberOfVectors*=(Integer)iter.next();  
        }
        
        Variable [] variables=new Variable[numberOfVariables];
        
        for(i=0; i<numberOfVectors; i++){
            
            for(int y=0; y<numberOfVariables; y++){
                variables[y]=new Real();
                variables[y].setValue((Double)((ArrayList)var.get(y)).get(index[y]));   // ??? kako lahko spremeni solution list
                
                if(y<1){
                    if(index[y]==vectorLengths.get(y)-1){
                        index[y]=0;
                        pogoj[y]=true;
                    }
                    else{
                        index[y]++;
                        pogoj[y]=false;
                    }   
                }
                
                if(y>0){
                    if(index[y]==vectorLengths.get(y)-1 && pogoj[y-1]){
                        index[y]=-1;
                        pogoj[y]=true;
                    }    
                                    
                    if(pogoj[y-1]){
                        index[y]++;
                        pogoj[y-1]=false;
                    }
                }
                
                //System.out.printf("X[%d]= % f ",y, variables[y].getValue());
            }
            
            //System.out.println();            
         
            newSolution = new Solution(problem_, variables);
            problem_.evaluate(newSolution);
            problem_.evaluateConstraints(newSolution);
            allSolutions.addToAll(newSolution); //Add solution to all-solution list
            ndl.add(newSolution);               //Add solution to non-dominated solution list 
        }
       
        allSolutions.printVariablesAndObjectivesToFile(logFileName);
        
        return ndl;
   }//execute
}//RandomSearch
