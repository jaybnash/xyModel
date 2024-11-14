import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

public class NModel {
    private double T;
    private int L;
    private int dim;
    private double[] spins; 
    private int totalSize;

    public NModel(double T, int L, int dim) {
        this.T = T;
        this.L = L;
        this.dim = dim;
        this.totalSize = (int) Math.pow(L, dim);
        this.spins = new double[totalSize];
        for (int i = 0; i < totalSize; i++) {
            spins[i] = Math.random() * 2 * Math.PI;
        }
    }

    private void reset() {
        for (int i = 0; i < totalSize; i++) {
            spins[i] = Math.random() * 2 * Math.PI;
        }
    }

    public void step() {
        double deltaTheta = 2 * Math.PI;
        int stepsPerCall = 10000; 

        for (int step = 0; step < stepsPerCall; step++) {
            int[] coords = getRandomCoords();
            int index = coordToIndex(coords);
            double thetaOld = spins[index];
            double deltaThetaRandom = (Math.random() - 0.5) * 2 * deltaTheta;
            double thetaNew = thetaOld + deltaThetaRandom;
            thetaNew = thetaNew % (2 * Math.PI);
            if (thetaNew < 0) thetaNew += 2 * Math.PI;

            double eDiff = deltaU(coords, thetaNew);

            if ((eDiff <= 0) || (Math.random() < Math.exp(-eDiff / this.T))) {
                spins[index] = thetaNew;
            }
        }
    }

    private double deltaU(int[] coords, double thetaNew) {
        double thetaOld = spins[coordToIndex(coords)];
        List<int[]> neighbors = getNeighbors(coords);
        double deltaU = 0.0;

        for (int[] neighborCoords : neighbors) {
            double neighborTheta = spins[coordToIndex(neighborCoords)];
            deltaU -= Math.cos(thetaNew - neighborTheta);
            deltaU += Math.cos(thetaOld - neighborTheta);
        }

        return deltaU;
    }

    private double magnetization() {
        double sumX = 0.0;
        double sumY = 0.0;
        for (double theta : spins) {
            sumX += Math.cos(theta);
            sumY += Math.sin(theta);
        }
        return Math.sqrt(sumX * sumX + sumY * sumY) / totalSize;
    }

    private double energy() {
        double energy = 0.0;
        for (int i = 0; i < totalSize; i++) {
            int[] coords = indexToCoord(i);
            List<int[]> neighbors = getNeighbors(coords);
            for (int[] neighborCoords : neighbors) {
                int neighborIndex = coordToIndex(neighborCoords);
                energy += -Math.cos(spins[i] - spins[neighborIndex]);
            }
        }
        // Each pair counted twice, so divide by 2
        return energy / (2 * totalSize);
    }

    private int countVortices() {
        if (dim != 2) {
            return 0;
        }
        int vortexCount = 0;
        // Assuming L >= 2 for vortices to exist
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                int right = (i + 1) % L;
                int bottom = (j + 1) % L;

                double theta0 = spins[coordToIndex(new int[]{i, j})];
                double theta1 = spins[coordToIndex(new int[]{right, j})];
                double theta2 = spins[coordToIndex(new int[]{right, bottom})];
                double theta3 = spins[coordToIndex(new int[]{i, bottom})];

                double dTheta = wrapAngle(theta1 - theta0) + wrapAngle(theta2 - theta1)
                              + wrapAngle(theta3 - theta2) + wrapAngle(theta0 - theta3);

                if (Math.abs(dTheta) >= 2 * Math.PI) {
                    vortexCount++;
                }
            }
        }
        return vortexCount;
    }

    private double wrapAngle(double angle) {
        angle = angle % (2 * Math.PI);
        if (angle < -Math.PI) angle += 2 * Math.PI;
        if (angle > Math.PI) angle -= 2 * Math.PI;
        return angle;
    }

    public Object[] simulate_energy_mag_sus(int L, float T_min, float T_max, int T_steps,
                                            int equilibration_steps, int measurement_steps, int dim, int n) {
        float dT = (T_max - T_min) / T_steps;
        float currentT = T_min;
        double[] energies = new double[T_steps];
        double[] magnetizations = new double[T_steps];
        double[] susceptibilities = new double[T_steps];
        int[] vortexCounts = new int[T_steps];
        
        for (int t = 0; t < T_steps; t++) {
            this.reset();
            this.T = currentT;
            // Equilibration
            for (int i = 0; i < equilibration_steps; i++) {
                this.step();
            }
            double E = 0.0;
            double E2 = 0.0;
            double M = 0.0;
            double M2 = 0.0;
            int vortexCount = 0;
            // Measurements
            for (int i = 0; i < measurement_steps; i++) {
                this.step();
                double energy = this.energy();
                double magnetization = this.magnetization();
                vortexCount += this.countVortices();
                E += energy;
                E2 += energy * energy;
                M += magnetization;
                M2 += magnetization * magnetization;
            }
            // Averaging
            E /= measurement_steps;
            E2 /= measurement_steps;
            M /= measurement_steps;
            M2 /= measurement_steps;
            energies[t] = E;
            magnetizations[t] = M;
            susceptibilities[t] = (this.dim == 2 ? (this.L * this.L) : Math.pow(this.L, this.dim)) * (1 / this.T) * (M2 - (M * M));
            vortexCounts[t] = (dim == 2) ? (vortexCount / measurement_steps) : 0;
            currentT += dT;
        }
        return new Object[]{energies, magnetizations, susceptibilities, vortexCounts};
    }

    public static void main(String[] args) {
        int[] latticeSizes = {10, 25, 50, 100};  // Different lattice sizes
        float T_min = 0.1f;
        float T_max = 2.5f;
        int T_steps = 25;
        int equilibration_steps = 100;
        int measurement_steps = 100;
        int n = 1; // Number of runs
        int dim = 5; // Number of dimensions

        String url = "jdbc:sqlite:./data/n_model.db";

        // Connect to the database
        try (Connection conn = DriverManager.getConnection(url)) {
            // Create table if it doesn't exist
            String createTableSQL = "CREATE TABLE IF NOT EXISTS SimulationResults ("
                    + "T REAL, "
                    + "L INTEGER, "
                    + "dim INTEGER, "
                    + "energy REAL, "
                    + "magnetization REAL, "
                    + "susceptibility REAL, "
                    + "vortex_count INTEGER)";
            conn.createStatement().execute(createTableSQL);

            String insertSQL = "INSERT INTO SimulationResults (T, L, dim, energy, magnetization, susceptibility, vortex_count) VALUES (?, ?, ?, ?, ?, ?, ?)";
            try (PreparedStatement pstmt = conn.prepareStatement(insertSQL)) {
                for (int L : latticeSizes) {
                    for (int run = 0; run < n; run++) {
                        NModel model = new NModel(T_min, L, dim);
                        Object[] results = model.simulate_energy_mag_sus(L, T_min, T_max, T_steps, equilibration_steps, measurement_steps, dim, n);
                        double[] energies = (double[]) results[0];
                        double[] magnetizations = (double[]) results[1];
                        double[] susceptibilities = (double[]) results[2];
                        int[] vortexCounts = (int[]) results[3];

                        // Insert each result set into the database
                        float currentT = T_min;
                        float dT = (T_max - T_min) / T_steps;
                        for (int i = 0; i < T_steps; i++) {
                            pstmt.setFloat(1, currentT);
                            pstmt.setInt(2, L);
                            pstmt.setInt(3, dim);
                            pstmt.setDouble(4, energies[i]);
                            pstmt.setDouble(5, magnetizations[i]);
                            pstmt.setDouble(6, susceptibilities[i]);
                            pstmt.setInt(7, vortexCounts[i]);
                            pstmt.executeUpdate();
                            currentT += dT;
                        }
                        System.out.println("Run " + (run + 1) + " completed for L=" + L + " in " + dim + "D");
                    }
                }
            }
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
    }

    // Helper method to convert n-dimensional coordinates to 1D index
    private int coordToIndex(int[] coords) {
        int index = 0;
        for (int i = 0; i < dim; i++) {
            index = index * L + coords[i];
        }
        return index;
    }

    // Helper method to convert 1D index to n-dimensional coordinates
    private int[] indexToCoord(int index) {
        int[] coords = new int[dim];
        for (int i = dim - 1; i >= 0; i--) {
            coords[i] = index % L;
            index /= L;
        }
        return coords;
    }

    // Helper method to get random coordinates in n dimensions
    private int[] getRandomCoords() {
        int[] coords = new int[dim];
        for (int i = 0; i < dim; i++) {
            coords[i] = (int) (Math.random() * L);
        }
        return coords;
    }

    // Helper method to get all neighbors of a given coordinate in n dimensions with periodic boundary conditions
    private List<int[]> getNeighbors(int[] coords) {
        List<int[]> neighbors = new ArrayList<>();
        for (int d = 0; d < dim; d++) {
            int[] neighborPos = coords.clone();
            neighborPos[d] = (coords[d] + 1) % L; // Positive direction
            neighbors.add(neighborPos);

            int[] neighborNeg = coords.clone();
            neighborNeg[d] = (coords[d] - 1 + L) % L; // Negative direction
            neighbors.add(neighborNeg);
        }
        return neighbors;
    }

    public double[] getSpins() {
        // Returns a copy to prevent external modification
        return spins.clone();
    }

    public double getT() {
        return T;
    }

    public void setT(double T) {
        this.T = T;
    }
}
