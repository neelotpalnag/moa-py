/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jmetal.metaheuristics.demoNN;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import jmetal.core.Algorithm;
import jmetal.core.Model;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.Kursawe;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;


/**
 *
 * @author Rok
 */

/**
 * Class for configuring and running the demo algorithm
 */
public class demoNN_main {
    public static Logger      logger_ ;      // Logger object
    public static FileHandler fileHandler_ ; // FileHandler object
    
    /**
   * @param args Command line arguments.
   * @throws JMException 
   * @throws IOException 
   * @throws SecurityException 
   * Usage: three choices
   *      - jmetal.metaheuristics.demo.demo_main
   *      - jmetal.metaheuristics.demo.demo_main problemName
   *      - jmetal.metaheuristics.demo.demo_mainproblemName paretoFrontFile
   */
    
   public static void main(String [] args) throws JMException, SecurityException, IOException, ClassNotFoundException {
       
        Problem   problem   ;         // The problem to solve
        Algorithm algorithm ;         // The algorithm to use
        Operator  selection ;
        Operator  crossover ;
        HashMap   parameters;         // Operator parameters
    
        QualityIndicator indicators ; // Object to get quality indicators
        
        // Logger object and file to store log messages
        logger_      = Configuration.logger_ ;
        fileHandler_ = new FileHandler("Demo_main.log"); 
        logger_.addHandler(fileHandler_) ;
        
        indicators = null ;
        if (args.length == 1) {
            Object [] params = {"Real"};
            problem = (new ProblemFactory()).getProblem(args[0],params);
        } // if
        else if (args.length == 2) {
            Object [] params = {"Real"};
            problem = (new ProblemFactory()).getProblem(args[0],params);
            indicators = new QualityIndicator(problem, args[1]) ;
        } // if
        else { // Default problem
            problem = new Kursawe("Real", 3); 
            //problem = new Water("Real");
              //problem = new ZDT1("ArrayReal", 100);
              //problem = new ConstrEx("Real");
              //problem = new DTLZ1("Real");
              //problem = new OKA2("Real") ;
         } // else
        
        algorithm = new DEMONNGP(problem);
        
        // Algorithm parameters
        algorithm.setInputParameter("populationSize",100);
        algorithm.setInputParameter("maxIterations",250);
        algorithm.setInputParameter("environmentalSelection","s");
        
        // Crossover operator 
        parameters = new HashMap();
        parameters.put("CR", 0.5);
        parameters.put("F", 0.5);
        crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);   
        
        // Add the operators to the algorithm
        parameters = null;
        selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters) ;

        algorithm.addOperator("crossover",crossover);
        algorithm.addOperator("selection",selection);
        
        // Execute the Algorithm 
        long initTime = System.currentTimeMillis();
        SolutionSet population = algorithm.execute();
        long estimatedTime = System.currentTimeMillis() - initTime;
        
        // Result messages 
        logger_.info("Total execution time: "+estimatedTime + "ms");
        logger_.info("Variables values have been writen to file VAR");
        population.printVariablesToFile("VAR");    
        logger_.info("Objectives values have been writen to file FUN");
        population.printObjectivesToFile("FUN");
  
        if (indicators != null) {
          logger_.info("Quality indicators") ;
          logger_.info("Hypervolume: " + indicators.getHypervolume(population)) ;
          logger_.info("GD         : " + indicators.getGD(population)) ;
          logger_.info("IGD        : " + indicators.getIGD(population)) ;
          logger_.info("Spread     : " + indicators.getSpread(population)) ;
          logger_.info("Epsilon    : " + indicators.getEpsilon(population)) ;  
        } // if        
   }//main 
}//GDE3_main
