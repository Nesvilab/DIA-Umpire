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
package MSUmpire.PSMDataStructure;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import com.compomics.util.experiment.biology.PTM;
import java.io.Serializable;
import java.text.DecimalFormat;

/**
 * PTM class
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class ModificationInfo implements Serializable{
    private static final long serialVersionUID = 25425698745L;

    public PTM modification;
    public float massdiff;
    public float mass;
    public String site;
    
    public String GetKey(){
        if (modification != null) {
            String tmp = new DecimalFormat("#.###").format(massdiff);
            return "[" + tmp + "(" + site + ")]";
        }
        return "";
    }
}
