/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jmetal.metaheuristics.exhaustiveSearch;

import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
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
 * Class for configuring and running the ExhaustiveSearch algorithm
 */
public class ExhaustiveSearch_main {
    public static Logger      logger_ ;      // Logger object
    public static FileHandler fileHandler_ ; // FileHandler object
    
    /**
   * @param args Command line arguments.
   * @throws JMException
   * @throws IOException
   * @throws SecurityException
   * Usage: three options
   *      - jmetal.metaheuristics.exhaustiveSearch.RandomSearch_main
   *      - jmetal.metaheuristics.exhaustiveSearch.RandomSearch_main problemName
   */
    
    public static void main(String [] args) throws
                                  JMException, SecurityException, IOException, ClassNotFoundException {
        
        Problem   problem   ;         // The problem to solve
        Algorithm algorithm ;         // The algorithm to use
        Operator  crossover ;         // Crossover operator
        Operator  mutation  ;         // Mutation operator
        Operator  selection ;         // Selection operator   
        
        QualityIndicator indicators ; // Object to get quality indicators
        
        // Logger object and file to store log messages
        logger_      = Configuration.logger_ ;
        fileHandler_ = new FileHandler("ExhaustiveSearch_main.log"); 
        logger_.addHandler(fileHandler_);
        
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
          problem = new Kursawe("Real", 2);     //The second parameter specifies the vector length (x1, x2,...xN)
          //problem = new Water("Real");
          //problem = new ZDT1("ArrayReal", 1000);
          //problem = new ZDT4("BinaryReal");
          //problem = new WFG1("Real");
          //problem = new DTLZ1("Real");
          //problem = new OKA2("Real") ;
        } // else
        
        algorithm = new ExhaustiveSearch(problem);
        
        // Algorithm parameters
        algorithm.setInputParameter("step", 1.0);
   
        // Execute the Algorithm
        long initTime = System.currentTimeMillis();
        SolutionSet population = algorithm.execute();
        long estimatedTime = System.currentTimeMillis() - initTime;

        // Result messages 
        logger_.info("Total execution time: "+estimatedTime + "ms");
        logger_.info("Objectives values have been writen to file FUN");
        population.printObjectivesToFile("FUN");
        logger_.info("Variables values have been writen to file VAR");
        population.printVariablesToFile("VAR");      
    
        if (indicators != null) {
          logger_.info("Quality indicators") ;
          logger_.info("Hypervolume: " + indicators.getHypervolume(population)) ;
          logger_.info("GD         : " + indicators.getGD(population)) ;
          logger_.info("IGD        : " + indicators.getIGD(population)) ;
          logger_.info("Spread     : " + indicators.getSpread(population)) ;
          logger_.info("Epsilon    : " + indicators.getEpsilon(population)) ;
        } // if                   
    } 
}
