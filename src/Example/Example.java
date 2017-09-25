/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Example;
import Algorithm.Equations;
import Jama.Matrix;

/**
 *
 * @author Wai Pai Lee
 * 
 * Example Equations:
 * cos(x) + cos(y) + cos(z) = 3/5
 * cos(3x) + cos(3y) + cos(3z) = 0
 * cos(5x) + cos(5y) + cos(5z) = 0
 */
public class Example implements Equations {
    final private int EQNS_LENGTH = 3;
    final private int VARIATE_NUM = 3;
    
    public Example() {}
    
    @Override
    public double[] f(double[] guess) {
        double [] outVal = new double[VARIATE_NUM];
        
        for(int i = 0; i < outVal.length; i++) {
            outVal[i] = valF(guess, i);
        }
        
        return outVal;
    }
    
    @Override
    public double valF(double[] guess, int eqnNum) {
        double result = 0;
        
        switch(eqnNum) {
            case 0:
                result = Math.cos(guess[0]) + Math.cos(guess[1]) + Math.cos(guess[2]) - (3.0/5.0);
                break;
                
            case 1:
                result = Math.cos(3*guess[0]) + Math.cos(3*guess[1]) + Math.cos(3*guess[2]);
                break;
                
            case 2:
                result = Math.cos(5*guess[0]) + Math.cos(5*guess[1]) + Math.cos(5*guess[2]);
                break;
                
            default:
                assert true;
        }
        
        return result * -1;
    }
    
    @Override
    public double derivative(double [] guess, int row, int col) {
        int location = row * VARIATE_NUM + col;
        double result = 1;
        
        switch(location) {
            case 1:
                result *= Math.sin(guess[0]);
                break;
                
            case 2:
                result *= Math.sin(guess[1]);
                break;
                
            case 3:
                result *= Math.sin(guess[2]);
                break;
            
            case 4:
                result *= (3 * Math.sin(3 * guess[0]));
                break;
                
            case 5:
                result *= (3 * Math.sin(3 * guess[1]));     
                break;
                
            case 6:
                result *= (3 * Math.sin(3 * guess[2]));
                break;
            
            case 7:
                result *= (5 * Math.sin(5 * guess[0]));
                break;
                
            case 8:
                result *= (5 * Math.sin(5 * guess[1]));
                break;
                
            case 9:
                result *= (5 * Math.sin(5 * guess[2]));
                break;
                
            default:
                assert true;
                
        }
        
        return result;
    }
    
    @Override
    public Matrix jacobian(double [] guess) {
        Matrix matrixJ = new Matrix(EQNS_LENGTH, VARIATE_NUM);
        
        for(int r = 0; r < EQNS_LENGTH; r++) {
            for(int c = 1; c <= VARIATE_NUM; c++) {
                matrixJ.set(r, c-1, (derivative(guess, r, c)));
            }
        }
        
        return matrixJ;
    }
    
    @Override
    public Matrix jacobianT(double [] guess) {
        Matrix matrixJ = new Matrix(EQNS_LENGTH, VARIATE_NUM);
        
        for(int r = 0; r < EQNS_LENGTH; r++) {
            for(int c = 1; c <= VARIATE_NUM; c++) {
                matrixJ.set(r, c-1, (derivative(guess, c-1, r+1)));
            }
        }
        
        return matrixJ;
    }
    
    @Override
    public int getEqnNum() {
        return EQNS_LENGTH;
    }
            
    @Override
    public int getVariateNum() {
        return VARIATE_NUM;
    }
}
