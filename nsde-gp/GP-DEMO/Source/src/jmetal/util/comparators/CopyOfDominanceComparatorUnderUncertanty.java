//  DominanceComparator.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.util.comparators;

import jmetal.core.Solution;

import java.util.Comparator;

/**
 * This class implements a <code>Comparator</code> (a method for comparing
 * <code>Solution</code> objects) based on a constraint violation test + 
 * dominance checking, as in NSGA-II.
 */
public class CopyOfDominanceComparatorUnderUncertanty implements Comparator{
 
  /** 
   * stores a comparator for check the OverallConstraintComparator
   */
  private static final Comparator overallConstraintViolationComparator_ =
                              new OverallConstraintViolationComparator();
 /**
  * Compares two solutions.
  * @param object1 Object representing the first <code>Solution</code>.
  * @param object2 Object representing the second <code>Solution</code>.
  * @return -1, or 0, or 1 if solution1 dominates solution2, both are 
  * non-dominated, or solution1  is dominated by solution22, respectively.
  * -2 sugests that solution1 could be better and that it should be exactly evaluated
  * 2 sugests that solution2 could be better and that it should be exactly evaluated
  * -10 sugests that solution2 could be better but is already exactly calculated so calculate again solution1
  * 10 sugests that solution1 could be better but is already exactly calculated so calculate again solution2
  * #100 we know nothing - not implemented
  */
  public int compare(Object object1, Object object2) {
    if (object1==null)
      return 1;
    else if (object2 == null)
      return -1;
    
    Solution solution1 = (Solution)object1;
    Solution solution2 = (Solution)object2;
        

    int dominate1 ; // dominate1 indicates if some objective of solution1 
                    // dominates the same objective in solution2. dominate2
    int dominate2 ; // is the complementary of dominate1.

    dominate1 = 0 ; 
    dominate2 = 0 ;
    
    int flag; //stores the result of the comparison

    if (solution1.getOverallConstraintViolation()!= 
        solution2.getOverallConstraintViolation() &&
       (solution1.getOverallConstraintViolation() < 0) ||         
       (solution2.getOverallConstraintViolation() < 0)){            
      return (overallConstraintViolationComparator_.compare(solution1,solution2));
    }
                                                
    // Equal number of violated constraints. Applying a dominance Test then
    double value1, value2;
    double variance1, variance2;
    int numberOfObjectiveswithTheSameValue = 0;
    
    for (int i = 0; i < solution1.numberOfObjectives(); i++) {
      value1 = solution1.getObjective(i);
      variance1 = solution1.getstandardDeviance(i);
      value2 = solution2.getObjective(i);
      variance2 = solution2.getstandardDeviance(i);
      
      if (value1 + variance1 < value2 - variance2) {
        flag = -1;
      } else if (value1 - variance1 > value2 + variance2) {
        flag = 1;
      } else {
        flag = 0;
        
        if(value1 + variance1 == value2 - variance2)
        	numberOfObjectiveswithTheSameValue ++;
        else if(value1 - variance1 > value2 + variance2)
        	numberOfObjectiveswithTheSameValue ++;
      }
      
      if (flag == -1) {
        dominate1++;
      }
      
      if (flag == 1) {
        dominate2++;           
      }
    }
            
    
    if(dominate1 > 0 && dominate2 > 0)
    	return 0;		//uncomparable
    
    
    if(dominate1 > 0 && dominate2 == 0)
    {
  	  if(dominate1 + numberOfObjectiveswithTheSameValue == solution1.numberOfObjectives())
  	  	return -1;
  	  else
  	  {
  	  	if(solution1.isExactllyEvaluated() == false)
  	  		return -2;
  	  	else
  	  		return 10;
  	  }
    } 
  	  
  	
    if(dominate1 == 0 && dominate2 > 0)
    {
  	  if(dominate2 + numberOfObjectiveswithTheSameValue == solution2.numberOfObjectives())
  	  	return 1;
  	  else
  	  {
  	  	if(solution2.isExactllyEvaluated() == false)
  	  		return 2;
  	  	else
  	  		return -10;
  	  }
    } 
 
    else		//(dominate1 == 0 && dominate2 == 0)
    {
  		if(solution1.isExactllyEvaluated() == true && solution2.isExactllyEvaluated() == true)
  			return 0;		//rešitvi imata enake vrednosti kriterijev    		
  		else	if (solution1.isExactllyEvaluated() == true)//becouse of variance we dont know which solution is better, we calculate one and compare again
  			return 2; 
  		else	//becouse of variance we dont know which solution is better, we calculate one and compare again
  			return -2; 
  	}
    
  } // compare
} // DominanceComparator


