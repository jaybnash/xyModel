import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;

// XYModel class for simulating the XY model using the metropolis algorithm
// Based an ising model implimentation provided by Dr. Daniel V. Schroeder: https://physics.weber.edu/thermal/isingJava.html
// Note: 
// This code does not include a visualization, as it is intended to be independent of any GUI or visualization code
// Additionally, this code relies upon having access to a sqlite-jbdc-3.47.0.0.jar library in the classpath
// The accompanying file Visualize.java can be used to visualize the XY model, and does not need the sqlite library in the classpath to complie and run
// Jay Nash 11/15/2024

public class XYModel {
    // Variables, T is the current temperature, L is the length of a side of the lattice, spins is the main 2D array holding spin values (0 to 2pi)
    private double T;
    private int L;
    private double[][] spins;

    // Constructor method, initializes T and L variables and sets a random value for each spin in the lattice
    public XYModel(double T, int L) {
        this.T = T;
        this.L = L;
        this.spins = new double[L][L];
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                spins[i][j] = Math.random() * 2 * Math.PI;
            }
        }
    }

    // Method to reset the lattice spins to random values, used instead of creating a new object each time
    private void reset() {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                spins[i][j] = Math.random() * 2 * Math.PI;
            }
        }
    }

    // Main step method for the XYModel, performs the metropolis algorithm to update spins
    public void step() {
        double deltaTheta = 2 * Math.PI; // Maximum change in angle for a spin (max value of possible spin values)
        
        for (int step=0; step<10000; step++) {         // update a randomly selected spin n times per step
            int i = (int) (Math.random() * this.L);    // row index selected randomly
            int j = (int) (Math.random() * this.L);    // column index selected randomly
            double thetaOld = this.spins[i][j];        // current spin value before change

            // to get the change in spin value we do the following:
            // generate a random number between 0.0 and 1.0 then subtract 0.5 to get a range of -0.5 to 0.5
            // multiply by 2 to get a range of -1 to 1, then multiply by deltaTheta to get a range of -2pi to 2pi
            // this means the spin can change in any direction up to a max change of 2pi
            double deltaThetaRandom = (Math.random() - 0.5) * 2 * deltaTheta;

            // new spin value is the old value plus the random change in spin
            double thetaNew = thetaOld + deltaThetaRandom;
            
            // new spin value must be between 0 and 2pi, modulo is used to adjust this
            thetaNew = thetaNew % (2 * Math.PI);

            // if we got a negative spin value, adjust it to be positive
            // this works because -0.5pi + 2pi = 1.5pi which is the same position on the unit circle (for example)
            if (thetaNew < 0) thetaNew += 2 * Math.PI;

            // calculate the change in energy for the new spin value
            double eDiff = deltaU(i, j, thetaNew);

            // the change in spin takes effect if energy decreases or if a random probability is met
            if ((eDiff <= 0) || (Math.random() < Math.exp(-eDiff / this.T))) this.spins[i][j] = thetaNew;
        }

    }

    // Method to calculate the difference in energy between the new spin value and the old spin value
    private double deltaU(int i, int j, double thetaNew) {
        double thetaOld = this.spins[i][j];                   // get old spin value
        double leftTheta, rightTheta, topTheta, bottomTheta;  // init variables for the 4 surrounding spins

        // get the surrounding spin values, wrapping around the array if at an edge
        if (i == 0) leftTheta = this.spins[this.L-1][j]; else leftTheta = this.spins[i-1][j];
        if (i == this.L-1) rightTheta = this.spins[0][j]; else rightTheta = this.spins[i+1][j];
        if (j == 0) topTheta = this.spins[i][this.L-1]; else topTheta = this.spins[i][j-1];
        if (j == this.L-1) bottomTheta = this.spins[i][0]; else bottomTheta = this.spins[i][j+1];

        // deltaU is calculated by:
        // first finding the new energy value by summing the cosines of the differences between the new spin value and the surrounding spins
        // second finding the old energy value by the same method
        // subtracting the old energy value by the new energy value to get the change in energy
        double deltaU = - (Math.cos(thetaNew - leftTheta) + Math.cos(thetaNew - rightTheta) +
                           Math.cos(thetaNew - topTheta) + Math.cos(thetaNew - bottomTheta))
                        + (Math.cos(thetaOld - leftTheta) + Math.cos(thetaOld - rightTheta) +
                           Math.cos(thetaOld - topTheta) + Math.cos(thetaOld - bottomTheta));
        return deltaU;
    }

    // Method to calculate the magnetization of the lattice
    private double magnetization() {
        double sumX = 0;  // init variable for sum of magnetic force in x direction
        double sumY = 0;  // init variable for sum of magnetic force in y direction

        // for each spin in the lattice, sum the x and y components of the magnetic force
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                sumX += Math.cos(-1 * spins[i][j]);
                sumY += Math.sin(-1 * spins[i][j]);
            }
        }

        // the total magnetic force is the square root of the sum of the squares of the x and y components
        // we then divide by L^2 to get the average magnitude of the magnetic force per spin
        // this method relys on the fact that we do not need the direction of the magnetic force, only the magnitude
        return Math.sqrt(sumX * sumX + sumY * sumY) / (L * L);
    }

    // Method to calculate the energy of the lattice
    private double energy() {
        double energy = 0;  // init variable for total energy

        // for each spin in the lattice, calculate the energy contribution from the surrounding spins
        // the formula used in this method follows directly from the method of deltaU
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                double leftTheta, rightTheta, topTheta, bottomTheta;
                if (i == 0) leftTheta = this.spins[this.L-1][j]; else leftTheta = this.spins[i-1][j];
                if (i == this.L-1) rightTheta = this.spins[0][j]; else rightTheta = this.spins[i+1][j];
                if (j == 0) topTheta = this.spins[i][this.L-1]; else topTheta = this.spins[i][j-1];
                if (j == this.L-1) bottomTheta = this.spins[i][0]; else bottomTheta = this.spins[i][j+1];
                energy += -Math.cos(spins[i][j] - leftTheta);
                energy += -Math.cos(spins[i][j] - rightTheta);
                energy += -Math.cos(spins[i][j] - topTheta);
                energy += -Math.cos(spins[i][j] - bottomTheta);
            }
        }

        // the total energy is divided by 2 to account for overcounting of spins
        // we then divide by L^2 to get the average energy per spin
        return (energy * 0.5) / (this.L * this.L);
    }

    // Method to count the number of vortices in the lattice
    private int countVortices() {
        int vortexCount = 0;  // init variable for total number of vortices

        // for each 2x2 square in the lattice, we do the following:
        // get the spin value of each spin in the square
        // find the change in angle between each spin and the next spin
        // if the change in angle is greater than or equal to 2pi, we have a vortex
        for (int i = 0; i < L-1; i++) {
            for (int j = 0; j < L-1; j++) {

                double theta0 = spins[i][j];
                double theta1 = spins[i+1][j];
                double theta2 = spins[i][j+1];
                double theta3 = spins[i+1][j+1];

                double dTheta = wrapAngle(theta1 - theta0) + wrapAngle(theta2 - theta1)
                              + wrapAngle(theta3 - theta2) + wrapAngle(theta0 - theta3);
                
                // take the absolute value of the change in angle, as dTheta can be between -4pi and 4pi
                if (Math.abs(dTheta) >= 2 * Math.PI) vortexCount++;
            }
        }
        return vortexCount;
    }

    // Method to wrap an angle to be between -pi and pi
    private double wrapAngle(double angle) {
        angle = angle % (2 * Math.PI);               // make sure angle is between 0 and 2pi, it should already be but lets be sure
        if (angle < -Math.PI) angle += 2 * Math.PI;  // if angle is less than -pi, add 2pi so it is between -pi and pi
        if (angle > Math.PI) angle -= 2 * Math.PI;   // if angle is greater than pi, subtract 2pi so it is between -pi and pi
        return angle;                                // return the wrapped angle
    }

    // Method to simulate a run on a range of temperatures for some lattice size
    // Returns an array of arrays containing the energy, magnetization, susceptibility, and vortex count for each temperature
    public Object[] simulate_energy_mag_sus(int L, float T_min, float T_max, int T_steps, int equilibration_steps, int measurement_steps) {
        // calculate the change in temperature per step, then set the initial temperature
        float dT = (T_max - T_min) / T_steps;
        float T = T_min;

        // init arrays to store the energy, magnetization, susceptibility, and vortex count for each temperature
        double[] energies = new double[T_steps];
        double[] magnetizations = new double[T_steps];
        double[] susceptibilities = new double[T_steps];
        int[] vortexCounts = new int[T_steps];

        // for each temperature do the following:
        // reset the lattice spins to random values
        // run the metropolis algorithm for a number of equilibration steps
        // then run the metropolis algorithm for a number of measurement steps
        // calculate the energy, magnetization, susceptibility, and vortex count for each measurement step
        // store the average values for each temperature for later analysis
        for (int t = 0; t < T_steps; t++) {
            this.reset();
            this.T = T;

            // equilibration steps
            for (int i = 0; i < equilibration_steps; i++) {
                this.step();  // step the simulation
            }

            // variables to store energy, energy^2 (unused), magnetization, magnetization^2, and vortex count
            double E = 0;
            double E2 = 0;
            double M = 0;
            double M2 = 0;
            int vortexCount = 0;

            // measurement steps
            for (int i = 0; i < measurement_steps; i++) {
                this.step();  // step the simulation
                double energy = this.energy();                // calculate energy for current state
                double magnetization = this.magnetization();  // calculate magnetization for current state
                vortexCount += this.countVortices();          // count vortices for current state
                E += energy;                                  // sum energy
                E2 += energy * energy;                        // sum energy^2
                M += magnetization;                           // sum magnetization
                M2 += magnetization * magnetization;          // sum magnetization^2
            }

            // store measurements
            M = M / measurement_steps;  // calculate average magnetization over all measurement steps
            E = E / measurement_steps;  // calculate average energy over all measurement steps
            energies[t] = E;            // store average energy
            magnetizations[t] = M;      // store average magnetization

            // calculate susceptibility
            // this is done via the formula: L^2 * (1 / T) * (M^2 - <M>^2)
            // this formula uses L^2 to scale the susceptibility to the lattice size
            // the 1 / T factor is used to scale the susceptibility to the temperature
            // the (M^2 - <M>^2) factor is used to calculate the variance of the magnetization
            susceptibilities[t] = (this.L * this.L) * (1 / this.T) * (M2 - (M * M));

            // store average vortex count over all measurement steps
            vortexCounts[t] = vortexCount / measurement_steps;

            // increment temperature for next run
            T += dT;
        }

        // return data to caller in an array of arrays, since only one object may be returned per method
        return new Object[]{energies, magnetizations, susceptibilities, vortexCounts};
    }

    // Main method to run the simulation for a range of lattice sizes and store collected results in database
    public static void main(String[] args) {
        // data is collected in the following way:
        // a number of lattice sizes are specified
        // temperatures to simulate are selected via a range and number of steps in that range
        // the number of equalibration steps and measurement steps are specified
        // the number of times to collect data for each lattice size is specified
        // a model is created, then simulate_energy_mag_sus is called for each lattice size
        // data is stored in a SQLite database for later analysis

        int[] latticeSizes = {10, 25, 50, 100, 200};    // Arbitrary list of lattice sizes to simulate
        float T_min = 0.1f;                             // Minimum temperature of runs
        float T_max = 2.5f;                             // Maximum temperature of runs
        int T_steps = 50;                               // Number of temperatures to simulate between min and max
        int equilibration_steps = 100;                  // Number of equilibration steps to run before measurements
        int measurement_steps = 100;                    // Number of measurement steps to run for each temperature
        int n = 10;                                     // Number of times to collect data for each lattice size

        String url = "jdbc:sqlite:./data/xy_model.db";  // Database path

        // Connect to sqlite database
        try (Connection conn = DriverManager.getConnection(url)) {
            // if the database does not exist, create it
            // the database stores the temp, lattice size, energy, magnetization, susceptibility, and vortex count for each run
            String createTableSQL = "CREATE TABLE IF NOT EXISTS SimulationResults ("
                    + "T REAL, "
                    + "L INTEGER, "
                    + "energy REAL, "
                    + "magnetization REAL, "
                    + "susceptibility REAL, "
                    + "vortex_count INTEGER)";
            conn.createStatement().execute(createTableSQL);

            // sql command prepared for inserting data into the database
            String insertSQL = "INSERT INTO SimulationResults (T, L, energy, magnetization, susceptibility, vortex_count) VALUES (?, ?, ?, ?, ?, ?)";

            // this try handles the insertion of data into the database via a prepared statement
            try (PreparedStatement pstmt = conn.prepareStatement(insertSQL)) {
                // for each lattice size, begin data collection
                for (int L : latticeSizes) {

                    // run the simulation n times
                    for (int run = 0; run < n; run++) {
                        // create a new model for each run
                        XYModel model = new XYModel(T_min, L);

                        // run simulation and get results
                        Object[] results = model.simulate_energy_mag_sus(L, T_min, T_max, T_steps, equilibration_steps, measurement_steps);
                        double[] energies = (double[]) results[0];
                        double[] magnetizations = (double[]) results[1];
                        double[] susceptibilities = (double[]) results[2];
                        int[] vortexCounts = (int[]) results[3];

                        // insert data into database, calculating the temperature for each run
                        // we can calculate the temperature now because we know how the temperature changes between runs and how the data is structured
                        // we insert the data for each temperature simulated
                        float T = T_min;
                        float dT = (T_max - T_min) / T_steps;
                        for (int i = 0; i < T_steps; i++) {
                            // setup prepared sql statement
                            pstmt.setFloat(1, T);
                            pstmt.setInt(2, L);
                            pstmt.setDouble(3, energies[i]);
                            pstmt.setDouble(4, magnetizations[i]);
                            pstmt.setDouble(5, susceptibilities[i]);
                            pstmt.setInt(6, vortexCounts[i]);

                            // execute prepared sql statement
                            pstmt.executeUpdate();

                            // increment temperature for next data point
                            T += dT;
                        }

                        // print completion message for each run
                        System.out.println("Run " + (run + 1) + " completed for L=" + L);
                    }
                }
            }
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
    }

    // Getters and setters for the XYModel class, spins should never be set so there is no setter for it
    public double[][] getSpins() {
        return spins;
    }

    public double getT() {
        return T;
    }

    public void setT(double T) {
        this.T = T;
    }
}
