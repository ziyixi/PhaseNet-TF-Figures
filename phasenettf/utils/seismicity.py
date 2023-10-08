import pandas as pd
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from pyproj import Proj
from typing import Callable
import warnings
from sklearn.exceptions import DataConversionWarning
from pandas.errors import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.filterwarnings(action="ignore", category=UserWarning)


def fit_slab_upper_interface_from_seismicity(
    data: pd.DataFrame,
    lon_center: float = -178,
    lat_center: float = -19,
    threshold=0.05,
) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """
    Fit a slab upper interface from seismicity data, which contains the following columns: longitude, latitude, depth
    """
    data_all = data.copy()
    # all longitude between 0 and 360
    seismicity_data = data_all[data_all.depth > 70]
    seismicity_data_above_150 = data_all[data_all.depth <= 150]
    print(
        len(seismicity_data),
        len(data_all),
        len(data_all[data_all.depth > 450]),
        len(seismicity_data_above_150),
    )

    seismicity_data.loc[seismicity_data.longitude < 0, "longitude"] += 360
    seismicity_data_above_150.loc[
        seismicity_data_above_150.longitude < 0, "longitude"
    ] += 360

    proj = Proj(f"+proj=sterea +lon_0={lon_center} +lat_0={lat_center} +units=km")
    # convert lonitude and latitude to x and y
    seismicity_data.loc[:, "x"], seismicity_data.loc[:, "y"] = proj(
        longitude=seismicity_data.longitude.values,
        latitude=seismicity_data.latitude.values,
    )
    seismicity_data_above_150.loc[:, "x"], seismicity_data_above_150.loc[:, "y"] = proj(
        longitude=seismicity_data_above_150.longitude.values,
        latitude=seismicity_data_above_150.latitude.values,
    )

    # fit a polynomial surface to the seismicity data, from x,y to z
    def inner_fit(df):
        model = make_pipeline(
            PolynomialFeatures(degree=4, include_bias=True), LinearRegression()
        )
        model.fit(X=df[["x", "y"]], y=df.depth.values)
        return model

    interface_step1 = inner_fit(seismicity_data)
    # calculate the std of the residuals for depth
    seismicity_data["depth_pred_step1"] = interface_step1.predict(
        X=seismicity_data[["x", "y"]]
    )
    seismicity_data["depth_residual_step1"] = (
        seismicity_data.depth - seismicity_data.depth_pred_step1
    )
    depth_residual_std = seismicity_data.depth_residual_step1.std()
    # remove outliers with abs depth difference>depth_residual_std*3
    seismicity_data = seismicity_data[
        np.abs(seismicity_data.depth_residual_step1) < depth_residual_std * 3
    ]
    # the second step, repeat the same process
    interface_step2 = inner_fit(seismicity_data)
    seismicity_data["depth_pred_step2"] = interface_step2.predict(
        X=seismicity_data[["x", "y"]]
    )
    seismicity_data["depth_residual_step2"] = (
        seismicity_data.depth - seismicity_data.depth_pred_step2
    )
    depth_residual_std = seismicity_data.depth_residual_step2.std()
    seismicity_data = seismicity_data[
        np.abs(seismicity_data.depth_residual_step2) < depth_residual_std * 3
    ]

    # do the same for the seismicity data above 150km
    interface_step1_above_150 = inner_fit(seismicity_data_above_150)
    seismicity_data_above_150["depth_pred_step1"] = interface_step1_above_150.predict(
        X=seismicity_data_above_150[["x", "y"]]
    )
    seismicity_data_above_150["depth_residual_step1"] = (
        seismicity_data_above_150.depth - seismicity_data_above_150.depth_pred_step1
    )
    depth_residual_std = seismicity_data_above_150.depth_residual_step1.std()
    seismicity_data_above_150 = seismicity_data_above_150[
        np.abs(seismicity_data_above_150.depth_residual_step1) < depth_residual_std * 3
    ]
    # the second step, repeat the same process
    interface_step2_above_150 = inner_fit(seismicity_data_above_150)
    seismicity_data_above_150["depth_pred_step2"] = interface_step2_above_150.predict(
        X=seismicity_data_above_150[["x", "y"]]
    )
    seismicity_data_above_150["depth_residual_step2"] = (
        seismicity_data_above_150.depth - seismicity_data_above_150.depth_pred_step2
    )
    depth_residual_std = seismicity_data_above_150.depth_residual_step2.std()
    seismicity_data_above_150 = seismicity_data_above_150[
        np.abs(seismicity_data_above_150.depth_residual_step2) < depth_residual_std * 3
    ]

    cur_ratio = 1.0
    step = 0
    length_seismicity_data = len(seismicity_data)
    while cur_ratio > threshold:
        # fit the surface again
        interface_cur = inner_fit(seismicity_data)
        depth_pred_cur = interface_cur.predict(X=seismicity_data[["x", "y"]])
        # remove all points with depth larger than depth_pred_cur
        seismicity_data = seismicity_data[seismicity_data.depth < depth_pred_cur]
        cur_ratio = len(seismicity_data) / length_seismicity_data
        print(f"step {step}, ratio={cur_ratio}, total={len(seismicity_data)}")
        step += 1

    # do the same for the seismicity data above 150km
    cur_ratio = 1.0
    step = 0
    length_seismicity_data = len(seismicity_data_above_150)
    while cur_ratio > threshold:
        # fit the surface again
        interface_cur_150 = inner_fit(seismicity_data_above_150)
        depth_pred_cur_150 = interface_cur_150.predict(
            X=seismicity_data_above_150[["x", "y"]]
        )
        # remove all points with depth larger than depth_pred_cur
        seismicity_data_above_150 = seismicity_data_above_150[
            seismicity_data_above_150.depth < depth_pred_cur_150
        ]
        cur_ratio = len(seismicity_data_above_150) / length_seismicity_data
        print(
            f"150km- step {step}, ratio={cur_ratio}, total={len(seismicity_data_above_150)}"
        )
        step += 1

    # return a function, it accepts np array of longitude and latitude, and return the depth
    # when the returned depth is smaller than 150km, use the interface_cur_150
    # when the returned depth is larger than 150km, use the interface_cur
    def interface(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
        x, y = proj(longitude=lon, latitude=lat)
        X_input = np.vstack([x, y]).T
        depth = interface_cur.predict(X=X_input)
        # depth_if_above_150 = interface_cur_150.predict(X=X_input)
        # # for all index with depth larger than 150km, use the same index value of depth_if_above_150
        # positions = np.where(depth <= 150)
        # depth[positions] = depth_if_above_150[positions]
        return depth

    return interface
