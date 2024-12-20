from pyology.enzymes import Enzyme

hexokinase = Enzyme(
    name="Hexokinase",
    k_cat=1.0,
    k_m={"glucose": 0.1},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["phosphoglucose_isomerase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

phosphoglucose_isomerase = Enzyme(
    name="Phosphoglucose Isomerase",
    k_cat=100,
    k_m={"glucose-6-phosphate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["phosphofructokinase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

phosphofructokinase = Enzyme(
    name="Phosphofructokinase",
    k_cat=100,
    k_m={"fructose-6-phosphate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["aldolase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

aldolase = Enzyme(
    name="Aldolase",
    k_cat=100,
    k_m={"fructose-1-6-bisphosphate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["triose_phosphate_isomerase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

triose_phosphate_isomerase = Enzyme(
    name="Triose Phosphate Isomerase",
    k_cat=100,
    k_m={"glyceraldehyde-3-phosphate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["glyceraldehyde-3-phosphate-dehydrogenase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

glyceraldehyde_3_phosphate_dehydrogenase = Enzyme(
    name="Glyceraldehyde 3-Phosphate Dehydrogenase",
    k_cat=100,
    k_m={"glyceraldehyde-3-phosphate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["phosphoglycerate_kinase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

phosphoglycerate_kinase = Enzyme(
    name="Phosphoglycerate Kinase",
    k_cat=100,
    k_m={"1-3-bisphosphoglycerate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["phosphoglycerate_mutase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

phosphoglycerate_mutase = Enzyme(
    name="Phosphoglycerate Mutase",
    k_cat=100,
    k_m={"3-phosphoglycerate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["enolase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

enolase = Enzyme(
    name="Enolase",
    k_cat=100,
    k_m={"2-phosphoglycerate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["pyruvate-kinase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

pyruvate_kinase = Enzyme(
    name="Pyruvate Kinase",
    k_cat=100,
    k_m={"phosphoenolpyruvate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["pyruvate-dehydrogenase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

pyruvate_dehydrogenase = Enzyme(
    name="Pyruvate Dehydrogenase",
    k_cat=100,
    k_m={"pyruvate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["pyruvate-carboxylase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

pyruvate_carboxylase = Enzyme(
    name="Pyruvate Carboxylase",
    k_cat=100,
    k_m={"pyruvate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["phosphoglycerate-mutate"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

phosphoglycerate_mutate = Enzyme(
    name="Phosphoglycerate Mutase",
    k_cat=100,
    k_m={"3-phosphoglycerate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["citrate-synthase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

citrate_synthase = Enzyme(
    name="Citrate Synthase",
    k_cat=100,
    k_m={"citrate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["aconitase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

aconitase = Enzyme(
    name="Aconitase",
    k_cat=100,
    k_m={"citrate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["isocitrate-dehydrogenase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

isocitrate_dehydrogenase = Enzyme(
    name="Isocitrate Dehydrogenase",
    k_cat=100,
    k_m={"isocitrate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["alpha-ketoglutarate-dehydrogenase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

alpha_ketoglutarate_dehydrogenase = Enzyme(
    name="Alpha-Ketoglutarate Dehydrogenase",
    k_cat=100,
    k_m={"alpha-ketoglutarate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["succinyl-coa-synthetase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

succinyl_coa_synthetase = Enzyme(
    name="Succinyl_CoA Synthetase",
    k_cat=100,
    k_m={"succinyl-coa": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["succinate-dehydrogenase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

succinate_dehydrogenase = Enzyme(
    name="Succinate Dehydrogenase",
    k_cat=100,
    k_m={"succinate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["fumarase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

fumarase = Enzyme(
    name="Fumarase",
    k_cat=100,
    k_m={"fumarate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["malate-dehydrogenase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)

malate_dehydrogenase = Enzyme(
    name="Malate Dehydrogenase",
    k_cat=100,
    k_m={"malate": 10},
    inhibitors={"adp": 1.0},
    activators={"atp": 1.0},
    active=True,
    downstream_enzymes=["oxaloacetate-synthetase"],
    hill_coefficients={"adp": 1.0, "atp": 1.0},
)
