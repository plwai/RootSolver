/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Algorithm;

import Jama.Matrix;

/**
 *
 * @author Wai Pai Lee
 */
public class NewtonMethod {
    final private int EQNS_LENGTH = 3;
    final private int VARIATE_NUM = 3;
    
    public NewtonMethod() {

    }
    
    public double[] f(double[] guess) {
        double [] outVal = new double[VARIATE_NUM];
        
        for(int i = 0; i < outVal.length; i++) {
            outVal[i] = valF(guess, i);
        }
        
        return outVal;
    }
    
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
    
    public double derivative(double [] guess, int row, int col) {
        int location = row * VARIATE_NUM + col;
        double result = -1;
        
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
    
    public Matrix jacobian(double [] guess) {
        Matrix matrixJ = new Matrix(EQNS_LENGTH, VARIATE_NUM);
        
        for(int r = 0; r < EQNS_LENGTH; r++) {
            for(int c = 1; c <= VARIATE_NUM; c++) {
                matrixJ.set(r, c-1, derivative(guess, r, c));
            }
        }
        
        return matrixJ;
    }
    
    public double[] solve(double [] guess) {
        double [] fx;
        double [] d = new double[VARIATE_NUM];
        double [] ans = new double[VARIATE_NUM];
        boolean loop = true;
        int itr = 0;
        int stop = 0;
        
        do {
            System.out.println(guess[0]);
            System.out.println(guess[1]);
            System.out.println(guess[2]);
            System.out.println();
            
            Matrix J = jacobian(guess);
            fx = f(guess);
            J.print(3, 10);
            //J.lu().getL().print(3, VARIATE_NUM);
            //J.lu().getU().print(3, VARIATE_NUM);
            (new Matrix(fx, VARIATE_NUM)).print(3, 10);
            if(J.lu().isNonsingular()) {
                d = J.lu().solve(new Matrix(fx, VARIATE_NUM)).getRowPackedCopy();
            }
            else {
                d[0] = 0.0001;
                d[1] = 0.0001;
                d[2] = 0.0001;
            }
            

            for(int i = 0; i < VARIATE_NUM; i++) {
                ans[i] = d[i] + guess[i];
            }
            
            for(int i = 0; i < VARIATE_NUM; i++) {
                if(Math.abs(d[i]) < 0.000001) {
                    loop = false;
                }
                else {
                    loop = true;
                    break;
                }
            }
            
            /*System.out.println(d[0]);
            System.out.println(d[1]);
            System.out.println(d[2]);
            System.out.println();*/
            
            
            guess = ans.clone();
            
            if(!loop) {
                stop++;
                loop = true;
            }
            else {
                stop = 0;
            }
            
            if(stop == 4) {
                break;
            }
            
            itr++;
        }while(loop);
                
        System.out.println(itr);
        
        return ans;
    }
    
    public static void main(String[] args) {
        NewtonMethod m = new NewtonMethod();
        
        double [] guess = {1,2,3};
        double [] ans;
        
        ans = m.solve(guess);
        
        System.out.println(ans[0]);
        System.out.println(ans[1]);
        System.out.println(ans[2]);
    }
}
