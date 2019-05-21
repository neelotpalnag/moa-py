//Author Miha Mlakar

package jmetal.encodings.solutionType;

import java.math.BigDecimal;
import java.text.DecimalFormat;

import jmetal.core.Problem;
import jmetal.core.SolutionType;
import jmetal.core.Variable;
import jmetal.encodings.variable.Real;
import jmetal.init.Params;


public class DiscreteSolutionType extends SolutionType {

	public DiscreteSolutionType(Problem problem) {
		super(problem);
		// TODO Auto-generated constructor stub
	}

	
	public Variable[] createVariables() throws ClassNotFoundException {
		Variable[] variables = new Variable[problem_.getNumberOfVariables()];
		Real real;
		double variable;		

		for (int var = 0; var < problem_.getNumberOfVariables(); var++)
		{
			real = new Real(problem_.getLowerLimit(var),problem_.getUpperLimit(var));   //MM
			variable = real.getValue();
			variable = getDiscretizedValue(variable, Params.getDiscretizationStep(var));	
			real.setValue(variable);
			variables[var] = real;
		}

		return variables ;
	}


	private double getDiscretizedValue(double variable, double discretizationStep) {

		double ldResult;
        int ldTempValue;
        int tempTemp =0;
        
		if (variable > 0)
            ldTempValue = (int)Math.floor((variable / discretizationStep) + 0.5);

        else
            ldTempValue = (int)Math.floor((variable / discretizationStep) - 0.5);

		//discretizationStep.
		BigDecimal b = new BigDecimal(discretizationStep);
		tempTemp = new Integer(ldTempValue);
        ldResult = (double)(ldTempValue * discretizationStep);
		
		return ldResult;
	}
	
	

}
