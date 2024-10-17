from pyology.common_enzymes import (
    aldolase,
    enolase,
    glyceraldehyde_3_phosphate_dehydrogenase,
    hexokinase,
    phosphofructokinase,
    phosphoglucose_isomerase,
    phosphoglycerate_kinase,
    phosphoglycerate_mutase,
    pyruvate_kinase,
    triose_phosphate_isomerase,
    phosphoglycerate_mutate,
)
from pyology.reaction import Reaction
from pyology.enzymes import Enzyme


class GlycolysisReactions:
    hexokinase = Reaction(
        name="Hexokinase",
        enzyme=Enzyme(name="Hexokinase", k_cat=1.0, k_m={"glucose": 0.1, "ATP": 0.1}),
        substrates={"glucose": 1, "ATP": 1},
        products={"glucose_6_phosphate": 1, "ADP": 1},
        reversible=False,
    )

    phosphoglucose_isomerase = Reaction(
        name="Phosphoglucose Isomerase",
        enzyme=phosphoglucose_isomerase,
        substrates={"glucose_6_phosphate": 1},
        products={"fructose_6_phosphate": 1},
        reversible=False,
    )

    phosphofructokinase = Reaction(
        name="Phosphofructokinase",
        enzyme=phosphofructokinase,
        substrates={"fructose_6_phosphate": 1, "ATP": 1},
        products={"fructose_1_6_bisphosphate": 1, "ADP": 1, "AMP": 1},
        reversible=False,
    )

    aldolase = Reaction(
        name="Aldolase",
        enzyme=aldolase,
        substrates={"fructose_1_6_bisphosphate": 1},
        products={"dihydroxyacetone_phosphate": 1, "glyceraldehyde_3_phosphate": 1},
        reversible=False,
    )

    triose_phosphate_isomerase = Reaction(
        name="Triose Phosphate Isomerase",
        enzyme=triose_phosphate_isomerase,
        substrates={"glyceraldehyde_3_phosphate": 1},
        products={"dihydroxyacetone_phosphate": 1},
        reversible=False,
    )

    phosphoglycerate_kinase = Reaction(
        name="Phosphoglycerate Kinase",
        enzyme=phosphoglycerate_kinase,
        substrates={"glyceraldehyde_3_phosphate": 1, "ADP": 1},
        products={"phosphoglycerate": 1, "ATP": 1},
        reversible=False,
    )

    phosphoglycerate_mutase = Reaction(
        name="Phosphoglycerate Mutase",
        enzyme=phosphoglycerate_mutate,
        substrates={"phosphoglycerate": 1},
        products={"dihydroxyacetone_phosphate": 1},
        reversible=False,
    )

    enolase = Reaction(
        name="Enolase",
        enzyme=enolase,
        substrates={"dihydroxyacetone_phosphate": 1},
        products={"phosphoenolpyruvate": 1},
        reversible=False,
    )

    pyruvate_kinase = Reaction(
        name="Pyruvate Kinase",
        enzyme=pyruvate_kinase,
        substrates={"phosphoenolpyruvate": 1, "ADP": 1},
        products={"pyruvate": 1, "ATP": 1},
        reversible=False,
    )

    glyceraldehyde_3_phosphate_dehydrogenase = Reaction(
        name="Glyceraldehyde 3-Phosphate Dehydrogenase",
        enzyme=glyceraldehyde_3_phosphate_dehydrogenase,
        substrates={"glyceraldehyde_3_phosphate": 1, "NAD": 1},
        products={"bisphosphoglycerate_1_3": 1, "NADH": 1},
        reversible=False,
    )
