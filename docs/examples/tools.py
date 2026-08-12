from xcompact3d_toolbox import Parameters
import traitlets
import os.path
import xarray as xr

VelocityDescription = {0: "Streamwise", 1: "Vertical", 2: "Spanwise"}


def get_datapath(case):
    return os.path.join("..", "incompact3d", "cases", case)

def prepare_dataset(prm, dataset):
    ds = invert_vertical_coordinate(dataset)
    ds = cut_to_test_section(prm, ds)
    clear_attrs_BC(ds)
    set_attrs(prm, ds)
    include_bed_position(prm, ds)
    return ds

def invert_vertical_coordinate(dataset):
    return dataset.assign_coords(y=dataset.y.data[::-1])

def cut_to_test_section(prm, dataset):
    ds = dataset.sel(x=slice(prm.xlx_pi, prm.xlx_pf))
    ds = ds.assign_coords(x=ds.x.data - prm.x0ramp)
    return ds

def clear_attrs_BC(dataset):
    for var in dataset.variables.keys():
        if "BC" in dataset[var].attrs:
            del dataset[var].attrs["BC"]

def include_bed_position(prm, dataset):
    y0ramp = 1.
    layer = prm.yly - prm.layer
    ramp = 1. + prm.declramp * dataset.x

    ramp[ramp < y0ramp] = y0ramp
    ramp[ramp > layer] = layer

    dataset['bed_position'] = ramp

    dataset.bed_position.attrs['name'] = 'Bed position'
    dataset.bed_position.attrs['long_name'] = '$x_{2r}$'

def set_attrs(prm, dataset):

    dataset.attrs[
        "title"
    ] = "Plunging condition for particle-laden flows over sloping bottoms: three-dimensional turbulence-resolving simulations"
    dataset.attrs["authors"] = "F. N. Schuch, E. Meiburg & J. H. Silvestrini"
    dataset.attrs["url"] = "https://github.com/fschuch/incompact3d_plunging_criterion"
    dataset.attrs["doi"] = "10.5281/zenodo.4044388"
    dataset.attrs["license"] = "Creative Commons Attribution 4.0 International"

    if "x" in dataset.coords:
        dataset.x.attrs = {"name": "Streamwise coordinate", "long_name": r"$x_1$"}
    if "y" in dataset.coords:
        dataset.y.attrs = {"name": "Vertical coordinate", "long_name": r"$x_2$"}
    if "z" in dataset.coords:
        dataset.z.attrs = {"name": "Spanwise coordinate", "long_name": r"$x_3$"}
    if "t" in dataset.coords:
        dataset.t.attrs = {"name": "Time", "long_name": r"$t$"}

    dataset["uset"] = xr.DataArray(
        prm.uset[0], attrs={"name": "Settling Velocity", "long_name": r"$u_s$"}
    )

    dataset["Ri0"] = xr.DataArray(
        prm.ri0, attrs={"name": "Initial Richardson Number", "long_name": r"$Ri_0$"}
    )

    dataset["Fr0"] = xr.DataArray(
        prm.fr0,
        attrs={"name": "Initial densimetric Froude Number", "long_name": r"$Fr_0$"},
    )

    dataset["Re"] = xr.DataArray(
        prm.re, attrs={"name": "Reynolds Number", "long_name": r"$Re$"}
    )

    dataset["S"] = xr.DataArray(
        prm.declramp, attrs={"name": "Bed Slop", "long_name": r"$S$"}
    )

def read_prm_to_dict(path):
    """
    This function reads the file BC-Plumes.prm from i3d and converts it into a Python dictionary
    
    Args:
        path: Path to folder that contains the file BC-Plumes.prm
    
    Returns:
        A dictionary containing flow parameters
    """
    prm = {}
    prm["filename"] = path + "/BC-Plumes.prm"
    prm["x0ramp"] = 0.0  # <- default value
    f = open(prm["filename"], "r")
    for line in f:
        line = line.partition("#")[0]  # Remove comments
        line = " ".join(line.split())  # Remove white spaces
        if line == "":  # Cycle if line is empty
            continue

        param = line.partition(" ")[0]
        value = line.partition(param + " ")[-1]
        if value[0] == "'" and value[-1] == "'":  # String
            value = value[1:-1]
        elif " " in value:  # Array
            if "." in value:
                value = [float(i) for i in value.split(" ")]  # Float
            else:
                value = [int(i) for i in value.split(" ")]  # int
        else:  # Not Array
            if "." in value:
                value = float(value)  # Float
            else:
                value = int(value)  # int
        prm[param] = value
    f.close()
    return prm


class PlumesParameters(Parameters):

    x0ramp = traitlets.Float(default_value=0.0).tag(
        group="Plumes", desc="horizontal position where ramp starts"
    )
    layer = traitlets.Float(default_value=1.0).tag(
        group="Plumes", desc="heigh of the flat plane downstream the ramp"
    )
    declramp = traitlets.Float(default_value=0.05).tag(
        group="Plumes", desc="slope declivity"
    )
    xlx_pi = traitlets.Float(default_value=0.0).tag(
        group="Plumes",
        desc="Physical domain start position (Dimensionless size in x-direction)",
    )
    xlx_pf = traitlets.Float(default_value=100.0).tag(
        group="Plumes",
        desc="Physical domain end position (Dimensionless size in x-direction)",
    )

    NameMappingForCompatibility = dict(nphi="numscalar", imodulo="ioutput")

    @property
    def ri0(self):
        if self.numscalar == 1:
            return self.ri[0]
        return self.ri

    @property
    def fr0(self):
        if self.numscalar == 1:
            return self.ri[0] ** (-0.5)
        return [ri0 ** (-0.5) for ri0 in self.ri]

    @property
    def datapath(self):
        return os.path.join(self.filename, "data")

    def load(self):

        dictionary = read_prm_to_dict(self.filename)

        # So values that are default for the plunging flow
        dictionary["nclyn"] = 1
        dictionary["filenamedigits"] = 1
        dictionary["ifilenameformat"] = "(I4.4)"

        for key in "ri uset cp".split():
            if not isinstance(dictionary[key], list):
                dictionary[key] = [dictionary[key]]

        # Boundary conditions are high priority in order to avoid bugs
        for bc in "nclx1 nclxn ncly1 nclyn nclz1 nclzn".split():
            if bc in dictionary:
                setattr(self, bc, dictionary[bc])

        for key, value in dictionary.items():
            Compatibilitykey = self.NameMappingForCompatibility.get(key, key)
            try:
                if self.trait_metadata(Compatibilitykey, "group") is not None:
                    setattr(self, Compatibilitykey, value)
            except:
                pass

