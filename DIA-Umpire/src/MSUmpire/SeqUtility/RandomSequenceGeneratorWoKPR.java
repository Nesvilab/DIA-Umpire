/*
 * This file is part of DIA-Umpire.
 *
 * DIA-Umpire is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * DIA-Umpire is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with DIA-Umpire. If not, see <https://www.gnu.org/licenses/>.
 *
 */
package MSUmpire.SeqUtility;

import java.util.Random;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;


/**
 *
 * @author Chih-Chiang Tsou
 */
public class RandomSequenceGeneratorWoKPR {

    private String RandomSeq;
    int RandomIdx = 0;
    /**
     * All possible characters
     */
    private final char[] CHARS = {
        'A', 'N', 'D', 'C', 'Q', 'E',
        'G', 'H', 'I', 'L', 'M', 'F',
        'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X'};

    public static RandomSequenceGeneratorWoKPR GetInstance(){
        if(randomSequenceGeneratorWoKPR==null){
            randomSequenceGeneratorWoKPR=new RandomSequenceGeneratorWoKPR();
        }
        return randomSequenceGeneratorWoKPR;
    }
    
    private static RandomSequenceGeneratorWoKPR randomSequenceGeneratorWoKPR=null;
    
    private RandomSequenceGeneratorWoKPR(){
        RandomSeq = generate(50);
    }
    
    transient ReadWriteLock lock = new ReentrantReadWriteLock();
   
    public char GetNext(){
        lock.writeLock().lock();
        try {
            if (RandomIdx == 50) {
                RandomIdx %= 50;
            }
            return RandomSeq.charAt(RandomIdx++);
        } finally {
            lock.writeLock().unlock();
        }
    }
    
    /**
     * Number of possible characters
     */
    private final int NUMBER_OF_CHARS = CHARS.length;
   
    
    /**
     * Random generator
     */
    private Random random = new Random();
    
            
    /**
     * Returns random sequence
     *
     * @param length Size of the sequence
     * @return Random sequence
     */
    private String generate(int length) {
        StringBuffer buffer = new StringBuffer();
        char randomChar;
        int randomInt;
        for (int i = 0; i < length; i++) {
            randomInt = random.nextInt(NUMBER_OF_CHARS);
            randomChar = CHARS[randomInt];
            buffer.append(randomChar);
        }
        return buffer.toString();
    }
}
