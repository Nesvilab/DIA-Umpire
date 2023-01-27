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

package MSUmpire.MathPackage;

/**
 *
 * @author Chih-Chiang Tsou
 */
public class MassDefect {

    public boolean InMassDefectRange(float mass, float d){
        //upper = 0.00052738*x + 0.066015 +0.1 
        //lower = 0.00042565*x + 0.00038210 -0.1

        double u = GetMassDefect(0.00052738d*mass + 0.066015d +d);
        double l = GetMassDefect(0.00042565d*mass + 0.00038210d -d);
        
        double defect=GetMassDefect(mass);
        if (u > l) {
            return (defect>=l && defect<=u);
        }
        return (defect>=l || defect<=u);
    }
    
    public double GetMassDefect(double mass){
        return mass-Math.floor(mass);
    }
    
}
