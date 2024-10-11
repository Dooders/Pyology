"""
Biological Background
Structure: 
    Double-membrane organelle with inner folds called cristae; 
    contains its own DNA.
Function:
    ATP Production: 
        Generates ATP through cellular respiration (glycolysis, Krebs cycle, 
        oxidative phosphorylation).
    Metabolic Integration: 
        Involved in apoptosis, calcium storage, and other metabolic pathways.
Modeling Considerations
    Metabolic Pathways:
        Glycolysis: 
        Occurs in the cytoplasm; glucose breakdown.
        Krebs Cycle: Occurs in the mitochondrial matrix.
        Electron Transport Chain (ETC): Located in the inner mitochondrial membrane.
    ATP Yield Calculations:
        Model the stoichiometry of ATP production from substrates.
"""

from organelle import Organelle


class Mitochondrion(Organelle):
    name = "Mitochondrion"
    """
    Represents the mitochondrion, a double-membrane organelle with inner folds
    called cristae that generates ATP through cellular respiration.

    Methods
    -------
    glycolysis : function
        Simulates the glycolysis of glucose.
    krebs_cycle : function
        Simulates the Krebs cycle processing of pyruvate.
    oxidative_phosphorylation : function
        Simulates the oxidative phosphorylation process.
    produce_atp : function
        Simulates the production of ATP from glucose.
    """

    def __init__(self) -> None:
        self.name = "Mitochondrion"
        self.atp = 0
        self.nadh = 0
        self.fadh2 = 0
        self.proton_gradient = 0
        self.proton_leak_rate = 0.1  # 10% of protons leak per cycle

    def electron_transport_chain(self) -> None:
        """
        Simulates the electron transport chain with more accurate proton pumping.
        """
        # Complex I: NADH to Ubiquinone
        if self.nadh > 0:
            protons_pumped = min(self.nadh * 4, 4)  # 4 H+ per NADH, max 4
            self.proton_gradient += protons_pumped
            self.nadh -= 1

        # Complex II: FADH2 to Ubiquinone (no protons pumped)
        if self.fadh2 > 0:
            self.fadh2 -= 1

        # Complex III: Ubiquinol to Cytochrome c
        protons_pumped = min(4, 2)  # 2 H+ per cycle, max 4
        self.proton_gradient += protons_pumped

        # Complex IV: Cytochrome c to Oxygen
        protons_pumped = min(4, 2)  # 2 H+ per 1/2 O2 reduced, max 4
        self.proton_gradient += protons_pumped

        print(f"Proton gradient after ETC: {self.proton_gradient}")

    def glycolysis(self, glucose_amount: int) -> int:
        """
        Simulates the glycolysis of glucose.
        """
        print(f"Glycolysis of {glucose_amount} units of glucose")
        # Simplify: Produce 2 ATP per glucose
        atp_produced = glucose_amount * 2
        self.atp += atp_produced
        pyruvate = glucose_amount * 2  # Each glucose yields 2 pyruvate molecules
        return pyruvate

    def krebs_cycle(self, pyruvate_amount: int) -> None:
        """
        Simulates the Krebs cycle processing of pyruvate.
        """
        print(f"Krebs cycle processing {pyruvate_amount} units of pyruvate")
        for _ in range(pyruvate_amount):
            self.atp += 1  # 1 ATP (or GTP) per pyruvate
            self.nadh += 4  # 3 NADH from Krebs cycle + 1 from pyruvate dehydrogenase
            self.fadh2 += 1  # 1 FADH2 per pyruvate

        print(
            f"After Krebs cycle: ATP: {self.atp}, NADH: {self.nadh}, FADH2: {self.fadh2}"
        )

    def oxidative_phosphorylation(self) -> None:
        """
        Simulates the oxidative phosphorylation process using the electron transport chain
        and ATP synthase, including proton leak.
        """
        print("Performing oxidative phosphorylation")

        # Run the electron transport chain
        self.electron_transport_chain()

        # Simulate proton leak
        leaked_protons = int(self.proton_gradient * self.proton_leak_rate)
        self.proton_gradient -= leaked_protons
        print(f"Protons leaked: {leaked_protons}")

        # Simulate ATP synthase using the proton gradient
        atp_produced = self.proton_gradient // 4  # Assume 4 protons needed per ATP
        self.atp += atp_produced
        self.proton_gradient %= 4  # Remaining protons after ATP synthesis

        print(f"ATP produced in this cycle: {atp_produced}")
        print(f"Remaining proton gradient: {self.proton_gradient}")

    def produce_atp(self, glucose_amount: int) -> int:
        pyruvate = self.glycolysis(glucose_amount)
        self.krebs_cycle(pyruvate)
        self.oxidative_phosphorylation()
        print(f"Total ATP produced: {self.atp}")
        return self.atp

    def function(self):
        """
        Implement the main function of the mitochondrion.
        """
        return "Produce energy through cellular respiration"


# Simulation code
if __name__ == "__main__":
    # Create a Mitochondrion instance
    mito = Mitochondrion()

    # Simulate ATP production with different amounts of glucose
    glucose_amounts = [1, 2, 5, 10]

    for glucose in glucose_amounts:
        print(f"\nSimulating ATP production with {glucose} glucose units:")
        mito.produce_atp(glucose)

        # Print the current state of the mitochondrion
        print(f"Current ATP: {mito.atp}")
        print(f"Remaining NADH: {mito.nadh}")
        print(f"Remaining FADH2: {mito.fadh2}")
        print(f"Remaining proton gradient: {mito.proton_gradient}")

        # Reset the mitochondrion for the next simulation
        mito.atp = 0
        mito.nadh = 0
        mito.fadh2 = 0
        mito.proton_gradient = 0

    print("\nSimulation complete.")
