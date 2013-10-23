/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.deciders;

/**
 *
 * @author rtahir
 */

    
    class DuplicatesException extends Exception {

    //Parameterless Constructor
    public DuplicatesException() {
        
        super ("DUPLICATES IN XML FILE");
    }

    //Constructor that accepts a message
    public DuplicatesException(String message) {
        super(message);
    }
}
    

