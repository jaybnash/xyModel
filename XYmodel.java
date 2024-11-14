import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;

public class XYModel {
    private double T;
    private int L;
    private double[][] spins;

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

    private void reset() {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                spins[i][j] = Math.random() * 2 * Math.PI;
            }
        }
    }

    public void step() {
        double deltaTheta = 2 * Math.PI;
        
        for (int step=0; step<10000; step++) {       // adjust number of steps as desired
            int i = (int) (Math.random() * this.L);    // choose a random row and column
            int j = (int) (Math.random() * this.L);
            double thetaOld = this.spins[i][j];
            double deltaThetaRandom = (Math.random() - 0.5) * 2 * deltaTheta; // random change in [-deltaTheta, deltaTheta]
            double thetaNew = thetaOld + deltaThetaRandom;
            thetaNew = thetaNew % (2 * Math.PI);
            if (thetaNew < 0) thetaNew += 2 * Math.PI;

            double eDiff = deltaU(i, j, thetaNew);

            if ((eDiff <= 0) || (Math.random() < Math.exp(-eDiff / this.T))) {    // Metropolis!
                this.spins[i][j] = thetaNew;
            }
        }

    }

    private double deltaU(int i, int j, double thetaNew) {
        double thetaOld = this.spins[i][j];
        double leftTheta, rightTheta, topTheta, bottomTheta;

        if (i == 0) leftTheta = this.spins[this.L-1][j]; else leftTheta = this.spins[i-1][j];
        if (i == this.L-1) rightTheta = this.spins[0][j]; else rightTheta = this.spins[i+1][j];
        if (j == 0) topTheta = this.spins[i][this.L-1]; else topTheta = this.spins[i][j-1];
        if (j == this.L-1) bottomTheta = this.spins[i][0]; else bottomTheta = this.spins[i][j+1];

        double deltaU = - (Math.cos(thetaNew - leftTheta) + Math.cos(thetaNew - rightTheta) +
                           Math.cos(thetaNew - topTheta) + Math.cos(thetaNew - bottomTheta))
                        + (Math.cos(thetaOld - leftTheta) + Math.cos(thetaOld - rightTheta) +
                           Math.cos(thetaOld - topTheta) + Math.cos(thetaOld - bottomTheta));
        return deltaU;
    }

    private double magnetization() {
        double sumX = 0;
        double sumY = 0;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                sumX += Math.cos(-1 * spins[i][j]);
                sumY += Math.sin(-1 * spins[i][j]);
            }
        }
        return Math.sqrt(sumX * sumX + sumY * sumY) / (L * L);
    }

    private double energy() {
        double energy = 0;
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
        return (energy * 0.5) / (this.L * this.L);
    }

    private int countVortices() {
        int vortexCount = 0;
        for (int i = 0; i < L-1; i++) {
            for (int j = 0; j < L-1; j++) {

                double theta0 = spins[i][j];
                double theta1 = spins[i+1][j];
                double theta2 = spins[i][j+1];
                double theta3 = spins[i+1][j+1];

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

    public Object[] simulate_energy_mag_sus(int L, float T_min, float T_max, int T_steps, int equilibration_steps, int measurement_steps) {
        float dT = (T_max - T_min) / T_steps;
        float T = T_min;
        double[] energies = new double[T_steps];
        double[] magnetizations = new double[T_steps];
        double[] susceptibilities = new double[T_steps];
        int[] vortexCounts = new int[T_steps];
        for (int t = 0; t < T_steps; t++) {
            this.reset();
            this.T = T;
            for (int i = 0; i < equilibration_steps; i++) {
                this.step();
            }
            double E = 0;
            double E2 = 0;
            double M = 0;
            double M2 = 0;
            int vortexCount = 0;
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
            M = M / measurement_steps;
            E = E / measurement_steps;
            energies[t] = E;
            magnetizations[t] = M;
            susceptibilities[t] = (this.L * this.L) * (1 / this.T) * (M2 - (M * M));
            vortexCounts[t] = vortexCount / measurement_steps;
            T += dT;
        }
        return new Object[]{energies, magnetizations, susceptibilities, vortexCounts};
    }

    public static void main(String[] args) {
        int[] latticeSizes = {10, 25, 50, 100, 200};  // Different lattice sizes
        float T_min = 0.1f;
        float T_max = 2.5f;
        int T_steps = 50;
        int equilibration_steps = 100;
        int measurement_steps = 100;
        int n = 10;

        String url = "jdbc:sqlite:./data/xy_model.db";

        // Connect to the database
        try (Connection conn = DriverManager.getConnection(url)) {
            // Create table if it doesn't exist
            String createTableSQL = "CREATE TABLE IF NOT EXISTS SimulationResults ("
                    + "T REAL, "
                    + "L INTEGER, "
                    + "energy REAL, "
                    + "magnetization REAL, "
                    + "susceptibility REAL, "
                    + "vortex_count INTEGER)";
            conn.createStatement().execute(createTableSQL);

            String insertSQL = "INSERT INTO SimulationResults (T, L, energy, magnetization, susceptibility, vortex_count) VALUES (?, ?, ?, ?, ?, ?)";
            try (PreparedStatement pstmt = conn.prepareStatement(insertSQL)) {
                for (int L : latticeSizes) {
                    for (int run = 0; run < n; run++) {
                        XYModel model = new XYModel(T_min, L);
                        Object[] results = model.simulate_energy_mag_sus(L, T_min, T_max, T_steps, equilibration_steps, measurement_steps);
                        double[] energies = (double[]) results[0];
                        double[] magnetizations = (double[]) results[1];
                        double[] susceptibilities = (double[]) results[2];
                        int[] vortexCounts = (int[]) results[3];

                        // Insert each result set into the database
                        float T = T_min;
                        float dT = (T_max - T_min) / T_steps;
                        for (int i = 0; i < T_steps; i++) {
                            pstmt.setFloat(1, T);
                            pstmt.setInt(2, L);
                            pstmt.setDouble(3, energies[i]);
                            pstmt.setDouble(4, magnetizations[i]);
                            pstmt.setDouble(5, susceptibilities[i]);
                            pstmt.setInt(6, vortexCounts[i]);
                            pstmt.executeUpdate();
                            T += dT;
                        }
                        System.out.println("Run " + (run + 1) + " completed for L=" + L);
                    }
                }
            }
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
    }

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
