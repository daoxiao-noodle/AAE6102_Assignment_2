# AAE6102_Assignment_2

Task 1, 4, and 5 are assisted by GenAI model, and the Chatroom Links are as follows:

Task 1: [https://chatgpt.com/share/680dd817-b934-8012-99cf-7513d533fd39](https://chatgpt.com/share/680dd817-b934-8012-99cf-7513d533fd39)

Task 4: [https://chatgpt.com/share/680de754-0294-8012-80c8-0deebe617d9e](https://chatgpt.com/share/680de754-0294-8012-80c8-0deebe617d9e)

Task 5:  [https://chatgpt.com/share/680df965-1168-8012-8196-c996258cd1ba](https://chatgpt.com/share/680df965-1168-8012-8196-c996258cd1ba)

# Task 1 – Comparative Analysis of Differential GNSS Positioning Techniques for Smartphone Navigation

## 1. Introduction

With the growing demand for high-precision navigation in smartphones, Global Navigation Satellite System (GNSS) technologies have evolved to offer better accuracy and reliability. Several advanced techniques are now being explored for smartphone integration, including Differential GNSS (DGNSS), Real-Time Kinematic (RTK), Precise Point Positioning (PPP), and PPP-RTK.

Each of these methods improves basic GNSS positioning in different ways:
- **Differential GNSS (DGNSS):** Enhances standard GNSS accuracy by using correction signals from ground-based reference stations or satellite-based augmentation systems, achieving meter-level improvements.
- **Real-Time Kinematic (RTK):** Provides centimeter-level positioning by using carrier-phase measurements and real-time corrections from nearby base stations.
- **Precise Point Positioning (PPP):** Achieves global high-accuracy positioning by using precise satellite orbit and clock corrections without relying on local base stations, although with longer convergence times.
- **PPP-RTK:** Combines the global accessibility of PPP with fast, centimeter-level convergence similar to RTK by integrating wide-area and regional correction services.

This essay first briefly introduces these techniques, then compares them systematically based on six critical dimensions: accuracy, infrastructure dependency, convergence and real-time usability, cost and practicality, hardware and complexity requirements, and environmental robustness.

## 2. Accuracy and Precision

- **DGNSS:** Enhances standalone GNSS accuracy from around 10 meters to approximately 1–5 meters, sufficient for general street navigation but inadequate for lane-level or high-precision needs.
- **RTK:** Provides exceptional centimeter-level accuracy (1–2 cm) by resolving carrier-phase ambiguities in real-time, making it suitable for applications such as lane positioning and augmented reality.
- **PPP:** Achieves decimeter to low-centimeter accuracy after convergence. Accuracy is better than DGNSS but generally less consistent than RTK under smartphone conditions.
- **PPP-RTK:** Combines global corrections with regional enhancements to reach centimeter-level precision with faster initialization, offering the best balance for mobile real-time navigation.

## 3. Infrastructure and Service Dependency

- **DGNSS:** Requires only basic augmentation signals, either broadcast via satellites or regional stations.
- **RTK:** Depends heavily on nearby base stations and real-time data communication through cellular networks.
- **PPP:** Needs only precise satellite correction data, typically available globally through satellite broadcast or internet distribution, without relying on local base stations.
- **PPP-RTK:** Depends on a network of regional and global correction sources, combining both satellite corrections and ground network assistance, requiring stable internet access.

## 4. Convergence Time and Real-Time Usability

- **DGNSS:** Provides immediate corrections, offering real-time usability without convergence delay—ideal for continuous smartphone navigation.
- **RTK:** Achieves convergence in seconds after receiving corrections, ensuring seamless, real-time centimeter-level tracking.
- **PPP:** Suffers from long convergence times, typically ranging from 5 to 30 minutes, which hinders real-time usability for smartphones.
- **PPP-RTK:** Shortens convergence to 1–5 minutes by resolving ambiguities quickly, making it practical for dynamic smartphone applications.

## 5. Cost and Practicality

- **DGNSS:** Free or included by default through public augmentation systems. No extra subscription or hardware cost.
- **RTK:** Involves substantial operational costs for correction service subscriptions and demands stable mobile data connectivity.
- **PPP:** Basic services may be free, but premium PPP services that offer faster convergence and higher accuracy may incur costs.
- **PPP-RTK:** Usually subscription-based, particularly for commercial-grade services (e.g., Trimble RTX, Qianxun SI). However, the cost is increasingly being bundled into smartphone service plans.

## 6. Hardware and Complexity Requirements

- **DGNSS:** Minimal additional hardware requirements. Current smartphones can easily leverage DGNSS corrections through standard single-frequency GNSS receivers.
- **RTK:** Requires dual-frequency (L1/L5) GNSS chips, advanced carrier-phase tracking algorithms, and robust antenna performance, which only high-end smartphones currently support.
- **PPP:** Requires good GNSS signal quality but is relatively hardware-friendly compared to RTK because it does not require ground-based reference station communication.
- **PPP-RTK:** Demands the most sophisticated integration, involving dual-frequency GNSS, continuous data reception, and rapid ambiguity resolution, challenging mid-range smartphone hardware.

## 7. Environmental Robustness

- **DGNSS:** Reasonably robust against mild obstructions but degrades in dense urban canyons due to multipath and signal blockage.
- **RTK:** Highly sensitive to signal interruptions, cycle slips, and multipath, requiring open-sky conditions for optimal performance.
- **PPP:** More tolerant to variable environments because of global corrections, but multipath and atmospheric conditions still impact convergence.
- **PPP-RTK:** Offers improved robustness over RTK by leveraging wide-area corrections but still benefits from open-sky views for best performance.

## 8. Conclusion

Each GNSS technique presents distinct advantages and limitations for smartphone navigation:

- **DGNSS** offers immediate and affordable improvements, suited for everyday applications.
- **RTK** provides the highest accuracy but faces barriers in cost, communication, and hardware complexity.
- **PPP** allows global high-accuracy navigation but is hindered by long convergence periods.
- **PPP-RTK** represents the best balance, offering fast, high-precision positioning with global coverage and is poised to become the mainstream solution for future smartphones.

As smartphone GNSS hardware continues to advance and correction service networks expand, PPP-RTK is expected to play a crucial role in enabling applications such as autonomous driving, pedestrian navigation, and augmented reality on handheld devices.

To further clarify the differences among the four GNSS techniques, the following table summarizes their key characteristics across critical performance dimensions.

## Summary Table

| Technology | Accuracy | Infrastructure Dependency | Convergence & Real-Time Usability | Cost & Practicality | Hardware & Complexity | Environmental Robustness |
|:---|:---|:---|:---|:---|:---|:---|
| **DGNSS** | 1–5 meters | Public augmentation systems | Immediate, real-time capable | Free or very low cost | Minimal hardware requirements | Good for general use, vulnerable to urban multipath |
| **RTK** | 1–2 cm | Local base station and network | Instant convergence if stable link | High cost for services | Requires dual-frequency GNSS and complex processing | Sensitive to signal blockage, needs open sky |
| **PPP** | Decimeter to low-centimeter after convergence | Global satellite corrections | Long convergence (5–30 min) | Basic services free; premium services paid | Moderate hardware demands | More robust than RTK; affected during convergence |
| **PPP-RTK** | Centimeter-level | Global + regional corrections | Fast convergence (1–5 min) | Moderate to high cost, often bundled | High-end GNSS receiver and processing needed | Better robustness than RTK, still favors open sky |

---

# Task 2 – GNSS in Urban Areas
To improve the GNSS positioning performance with the skymask provided, and test with the "Urban" data in Assignment 1. The corresponding implementation can be found in the `tracking_POS_WLS_TASK2.m` file.

As the skymask can identify satellite visibility blockage by providing satellite azimuth and elevation angles, the given skymask is visualized as follows, in which the obscured portion is shown in gray and the unobscured portion is shown in white:

<div align="center">
    <img src="/Figure/Task2_skymask.jpg" width="300">
</div>

In Assignment 1, the satellite for the urban dataset is decoded only 4 satellites, while at least 4 satellites are required for the GNSS positioning solution. Finding that problem, in this task, the code for assignment 1 is firstly corrected, and then decoded to 5 satellites as shown below:

<div align="center">
    <img src="/Figure/Task2_satellite_position.jpg" width="300">
</div>

Overlaying the decoded satellite positions with the sky mask, the satellite passes can be obtained as shown in the figure below, where the obscured satellites are shown in red and the unobscured satellites are shown in green:

<div align="center">
    <img src="/Figure/Task2_satellite_position_with_skymask.jpg" width="300">
</div>

As can be seen from the above figure, for the Urban dataset, only two satellites are unobstructed, while the other three satellites are all obstructed. As mentioned earlier, GNSS positioning solution requires at least 4 satellites, so it is not possible to eliminate these three satellites directly. Based on the above problem, the satellite weights are adjusted based on the sky mask in the weighted least squares(WLS) method, which means if the current satellite is occluded, it is multiplied by another weight on top of its original weights including **1.0, 0.8, 0.5, 0.3, 0.0**, while the unoccluded satellites retain full weight. The position results on geographic map is shown as follows:

<div align="center">
    <img src="/Figure/Task2_position_results_on_geographic_map.jpg" width="800">
</div>

And RMSE between estimated position and ground truth in the ENU coordinate system is shown as follows:

## Summary Table

| **Occluded Satellite Weight Coefficient** | **Unoccluded Satellite Weight** | **RMSE (meters)** | **Positioning Performance** | **Observations** |
|:---|:---|:---|:---|:---|
| **1.0** | 1.0 | 148.8942 | Poor | No distinction between occluded and unoccluded satellites; severe errors. |
| **0.8** | 1.0 | 120.6495 | Improved | Moderate reduction in occluded satellite impact; significantly better accuracy. |
| **0.5** | 1.0 | 109.7933 | Best | Optimal balance: suppressed noise from occluded satellites while maintaining good geometry. |
| **0.3** | 1.0 | 109.8837 | Slightly Worse | Over-suppression leads to no further RMSE improvement; slight geometry loss effect. |
| **0.0 (Excluded)** | 1.0 | 112.9278 | Degraded | Full exclusion worsens geometry (PDOP), causing RMSE increase despite clean data. |

## Analysis

### 1. Significant Impact of Occluded Satellites
Without distinguishing between occluded and unoccluded satellites (i.e., assigning a weight of 1.0 to all satellites), the RMSE reached the highest value of 148.8942 meters. This result indicates that treating all satellites equally, regardless of occlusion, introduces significant positioning errors. The degraded signals from occluded satellites severely contaminate the position estimation when no proper weighting is applied.

### 2. Effectiveness of Weight Adjustment
By reducing the weight of occluded satellites to 0.8, the RMSE notably improved to 120.6495 meters, demonstrating that moderately devaluing occluded signals can significantly enhance positioning accuracy. Further lowering the weight to 0.5 achieved the best result, reducing RMSE to 109.7933 meters, suggesting that a more aggressive but controlled down-weighting yields an optimal balance between usable information and noise suppression.

### 3. Over-suppression or Exclusion is Counterproductive
When the occluded satellite weights were further reduced to 0.3 or fully excluded (0.0), the RMSE did not continue to decrease. In fact, a slight degradation in performance was observed, with RMSE values slightly increasing compared to the best-case scenario. This phenomenon is attributed to the deterioration of satellite geometry (higher PDOP) caused by removing too many satellites, which outweighs the benefits of noise reduction.

### 4. Optimal Strategy
The optimal strategy is to apply a moderate weight reduction (around 0.5) to occluded satellites. This approach suppresses the adverse effects of signal degradation while preserving the valuable geometric information necessary for accurate positioning. It achieves the best compromise between minimizing noise influence and maintaining favorable satellite constellation geometry.

In GNSS-challenging environments such as urban canyons or semi-indoor spaces, it is not advisable to simply exclude all occluded satellites. Instead, it is recommended to adjust satellite weights based on occlusion indicators, such as sky mask visibility or carrier-to-noise ratio (CN₀). By retaining partially occluded satellites with appropriately reduced weights, positioning robustness can be significantly improved without sacrificing satellite coverage and geometry quality.

## Conclusion

When there are only a few available satellites with partial occlusion, maintaining full weights for unoccluded satellites and moderately reducing occluded satellite weights (around 0.5) significantly improves WLS positioning accuracy, outperforming the strategy of directly excluding occluded satellites.

---

# Task 3 – GPS RAIM

## 1. Introduction
Receiver Autonomous Integrity Monitoring (RAIM) is a critical technique used to detect and exclude faulty satellite measurements, ensuring the reliability of GNSS positioning without relying on external augmentation systems. In GNSS, signal anomalies may occur due to satellite failures, multipath, or signal blockage. RAIM enables receivers to autonomously monitor and remove erroneous measurements.

In this task, a Weighted Least Squares (WLS) + Weighted RAIM algorithm is implemented, with C/N₀-based weighting and iterative fault exclusion. The corresponding implementation can be found in the `tracking_POS_WLS_TASK3.m` file.

## 2. Methodology
### 2.1 Weighted Least Squares (WLS) Positioning

The pseudorange observation model:

$$
\mathbf{ρ_i} = ||\mathbf{x}_ {sv,i} - x_{usr}|| + c \cdot δt_{usr} + ε_i
$$

where:

- $ρ_i$ : Measured pseudorange to satellite ($i$)
- ${x}_{sv,i}$ : Satellite position
- ${x}_{usr}$ : User receiver position ($x$, $y$, $z$)
- $c$ : Speed of light
- $δ t_{usr}$ : Receiver clock bias
- $ε_i$ : Measurement noise

Linearized:

$$
\mathbf{y} = \mathbf{H} \cdot \mathbf{x} + \mathbf{v}
$$

where:

- $y$ : Residual vector
- $H$ : Design matrix
- $x$ : State vector [dx, dy, dz, dclock]
- $v$ : Measurement noise vector

The WLS solution is:

$$
\hat{\mathbf{x}} = (\mathbf{H}^T \mathbf{W} \mathbf{H})^{-1} \mathbf{H}^T \mathbf{W} \mathbf{y}
$$

where $W$ = $diag$ ( ${w}_ {1}$, ${w}_ {2}$, …, ${w}_ {n}$ ) is a diagonal weighting matrix.
Typically, ${w}_ {i}$  ∝ $C / {N}_ {0}$ or inversely proportional to expected measurement noise variance.

## 2.2 RAIM Fault Detection and Exclusion

The residual after WLS estimation:

$$
\mathbf{r} = \mathbf{y} - \mathbf{H} \hat{\mathbf{x}}
$$

### Test statistic:

$$
T = \mathbf{r}^T \mathbf{W} \left( \mathbf{I} - \mathbf{H} (\mathbf{H}^T \mathbf{W} \mathbf{H})^{-1} \mathbf{H}^T \mathbf{W} \right) \mathbf{r}
$$

Fault detection rule:
- If  $T$ > $Threshold$ , a fault is detected.
- Otherwise, the measurement set is considered fault-free.

### Threshold Setting

Given:
- Measurement noise standard deviation  σ = 3 m 
- Probability of missed detection  $P_{md}$ = $10^{-7}$

The RAIM detection threshold:

$$
\text{Threshold} = (5.33 · \sigma)^2
$$

Here, 5.33 is the number of standard deviations corresponding to $P_{md}$ = $10^{-7}$.

### Faulty Measurement Exclusion

If a fault is detected:
- Identify the satellite with the maximum residual.
- Exclude the corresponding measurement.
- Recalculate WLS solution.
- Repeat RAIM test until no fault is detected or fewer than 4 satellites remain.

## 2.3 3D Protection Level (PL) Computation

The 3D Protection Level quantifies the worst-case position error under specified integrity risk.

Given the WLS solution covariance:

$$
\mathbf{Q} = \left( \mathbf{H}^T \mathbf{W} \mathbf{H} \right)^{-1}
$$

Extract the position-related submatrix $\mathbf{Q}_{\text{pos}}$ ( top-left 3 × 3 block ).

The 3D Protection Level (PL) is computed as:

$$
{PL}_ {3D} = K_{md} · \sqrt{{λ}_{max}}  ·  σ
$$

Where:
- $K_{md}$ = 5.33 (for $P_{md}$ = $10^{-7}$ )
- ${λ}_ {max}$ is the largest eigenvalue of $\mathbf{Q}_{\text{pos}}$
- σ = 3 m

## 2.4 Stanford Chart Analysis

The Stanford Chart visualizes GNSS integrity performance:

- **X-axis**: Positioning Error (meters)
- **Y-axis**: Protection Level (meters)
- **Alarm Limit (AL)**: set to 50 meters

Classification:

| Condition           | Rule |
|:--------------------|:----|
| **Good Detection**   | ( Error > AL) ∧ ( PL > AL ) |
| **Missed Detection** | ( Error > AL ) ∧ ( PL ≤ AL ) |
| **False Alarm**      | ( Error ≤ AL ) ∧ ( PL > AL ) |
| **Correct No Alarm** | ( Error ≤ AL ) ∧ ( PL ≤ AL ) |

Each point on the chart represents an epoch and its classification outcome.

## 3. Experimental Results and Analysis
### 3.1 Positioning Results
The spatial distribution of positioning results is visualized on the geographic map below, in which the blue points represent raw WLS (without RAIM) solutions, while the green points represent WLS solutions after RAIM-based exclusion, and the red point denotes the known reference ("best") position.

<div align="center">
    <img src="/Figure/Task3_position_results_on_geographic_map.jpg" width="800">
</div>

#### Findings and Analysis
- **Improvement in Certain Regions:**
Compared to the WLS-only results, RAIM-excluded positions show reduced scattering in some areas. The green points are noticeably tighter around the true location compared to the widespread blue points in open spaces, indicating partial effectiveness of the RAIM fault exclusion process.

- **Persistent Positioning Errors:**
Despite some improvements, significant deviations remain even after RAIM processing. Many green points are still located hundreds of meters away from the ground truth, confirming that undetected faulty measurements continued to affect the positioning solution.

- **Localized Effectiveness:**
RAIM appears to perform better in certain areas (e.g., near open fields or less obstructed zones), but fails to maintain performance across the entire environment. This suggests that the method's fault detection sensitivity is highly environment-dependent and possibly related to satellite geometry and signal obstruction variations.

### 3.2 GNSS integrity monitoring performance using a Stanford Chart analysis

In this experiment, the Stanford Chart provides a graphical distribution of each epoch by plotting the actual Position Error versus the corresponding Protection Level (PL), as shown below, in which the green circles represent epochs categorized as Correct No Alarm (both Position Error and PL below the Alarm Limit), and the black crosses represent Missed Detections (Position Error above the Alarm Limit while PL remains underestimated):

<div align="center">
    <img src="/Figure/Task3_Stanford_Chart.jpg" width="800">
</div>

#### Findings and Analysis: 
- **Dominance of Missed Detection Points:**
The majority of the plotted points are black crosses, corresponding to missed detections. Most epochs with large positioning errors (up to 1200 meters) are accompanied by Protection Levels significantly below the Alarm Limit. This suggests that the Protection Level estimates are systematically underestimated and fail to reflect the actual navigation risk.

- **Clustering of Correct No Alarms:**
A small cluster of green circles (Correct No Alarm) appears close to the origin. These represent epochs where both the Position Error and the Protection Level remained below the Alarm Limit. However, these points are relatively few and mainly occur when the signal conditions are truly optimal.

- **Underperforming RAIM Detection Capability:**
The general distribution confirms that the RAIM implementation does not effectively adjust Protection Levels in response to varying signal conditions, thereby missing critical integrity threats.

Then based on this Stanford Chart categorization, the integrity events were quantitatively summarized as follows:

| Integrity Category     | Epochs Count |
|------------------------|--------------|
|  Correct No Alarm     | 128          |
|  Missed Detection     | 3205         |
|  False Alarm          | 0            |
|  Good Detection       | 0            |

#### Findings and Analysis: 
- **Severe Missed Detection Issue:**
Out of a total of 3333 epochs, only 128 epochs were correctly monitored without triggering any false alarms or missed detections. In contrast, 3205 epochs suffered from missed detection, accounting for approximately 96% of the total. This indicates that in the majority of epochs, the actual position error exceeded the defined Alarm Limit (50 meters) without corresponding growth in the Protection Level (PL), leading to undetected faults. Such a high missed detection rate suggests that the current RAIM implementation is unable to reliably detect measurement faults under the given conditions.

- **No Good Detections or False Alarms:**
The absence of both Good Detections and False Alarms implies that while the RAIM system did not incorrectly raise alarms when they were not needed, it also failed to recognize actual integrity threats. This behavior reflects a highly conservative or ineffective detection threshold, resulting in protection levels that do not adequately bound true positioning errors.


## 6. Conclusion and Discussion
- **Observable Improvements Over Baseline WLS:**
Compared with the original Weighted Least Squares (WLS) method, the RAIM-integrated solution demonstrates localized improvements in positioning accuracy. In multiple segments of the trajectory, the application of fault detection and exclusion led to the suppression of outlier estimates and partial recovery of the true path. This is evident in the geographic distribution, where the green points (RAIM-filtered) appear more clustered and less scattered than the original blue WLS results. In this sense, the residual-based exclusion strategy was partially successful in identifying and removing certain high-impact erroneous satellite measurements.

- **Performance Degradation Due to Excessive Satellite Exclusion:**
However, further inspection reveals that the RAIM algorithm frequently excludes up to four satellites from the original eight available. While this aggressive response to residual anomalies helps suppress extreme outliers and superficially improves the dispersion of positioning results, it drastically reduces the remaining measurement redundancy. With only four satellites retained—the minimum required for solving three-dimensional positioning with clock bias estimation—the system loses the degrees of freedom necessary for reliable fault detection and protection level adjustment. Consequently, the satellite geometry deteriorates severely, leading to poor Dilution of Precision (DOP) and ill-conditioned estimation matrices. This weakened geometric configuration compromises the positioning solution, resulting in systematic biases and significant position deviations even in the absence of major residual errors. The overall effect is confirmed by the extremely high missed detection rate observed in the Stanford Chart analysis, highlighting the critical trade-off between aggressive fault exclusion and the preservation of observability.

- **Inadequate Protection Level Boundaries:**
The Protection Level (PL), which should conservatively bound position error under worst-case assumptions, remains severely underestimated in most epochs. Despite actual errors exceeding 500 to 1000 meters in some cases, the computed PL often remains well below the 50-meter alarm threshold. This suggests that the current PL formulation lacks sensitivity to residual dispersion and measurement covariance variation, especially under degraded satellite geometry.

- **Limitations and Constraints in Implementation:**
It is worth noting that due to project time constraints, further refinement of the RAIM exclusion logic and PL computation could not be completed. The observed performance shortfalls are not necessarily indicative of the RAIM principle’s ineffectiveness, but rather of the current configuration's limitations—particularly the unrestricted exclusion strategy and static weighting model. More sophisticated mechanisms, such as capping the number of allowable exclusions or dynamically scaling residual thresholds based on geometry strength, could likely yield more balanced results.

- **Summary and Outlook:**
In summary, while the implemented RAIM method improves upon standard WLS in selected scenarios, its effectiveness is compromised by over-aggressive satellite exclusion and insufficient PL estimation under limited redundancy. Future enhancements should aim to maintain observability, prevent over-filtering, and more accurately reflect risk in protection bounds. These improvements are necessary to achieve the level of integrity assurance expected in safety-critical GNSS applications.

---

# Task 4 – Difficulties and challenges of Using LEO communication Satellites for GNSS Navigation

## 1. Introduction

Low Earth Orbit (LEO) satellites, typically operating between 160 km and 2,000 km above the Earth’s surface, are extensively utilized for communication services such as global internet coverage (e.g., Starlink, OneWeb). Recently, there has been growing interest in leveraging LEO constellations for navigation purposes, either as an enhancement or an alternative to traditional Medium Earth Orbit (MEO) GNSS systems like GPS, Galileo, and BeiDou. While LEO satellites offer promising benefits, including lower signal latency, higher received power, and more rapid updates, they also pose significant technical, operational, and regulatory challenges. This essay systematically discusses the major difficulties encountered in adapting LEO communication satellites for GNSS navigation.

## 2. High Orbital Dynamics and Coverage Limitations

**(a) Rapid Satellite Movement and Short Visibility**

LEO satellites move at approximately 7.8 km/s, completing an orbit in roughly 90 minutes. As a result, a single satellite remains visible to a ground user for only 10–20 minutes, necessitating frequent handovers. Unlike MEO GNSS satellites, which can be tracked for several hours, LEO-based navigation requires receivers to rapidly acquire, track, and switch among satellites, imposing high demands on tracking algorithms and receiver stability.

**(b) Need for Large Constellations**

Due to their low altitude, LEO satellites cover much smaller areas on the Earth's surface compared to MEO satellites. Providing continuous, global navigation services would thus require deploying hundreds or even thousands of satellites. For instance, Starlink plans to operate over 40,000 satellites to achieve global internet coverage; adapting similar scales for navigation purposes would substantially increase costs, operational complexity, and maintenance burdens.

## 3. Signal Processing Challenges

**(a) Severe Doppler Shift Effects**

The high velocity of LEO satellites results in substantial Doppler frequency shifts, often several kHz, much higher than those observed in traditional GNSS systems. Navigation receivers must be capable of continuous Doppler tracking and rapid frequency compensation, requiring more sophisticated hardware and increased computational resources.

**(b) Signals Not Optimized for Navigation**

LEO satellites primarily transmit broadband communication signals, which are not inherently suitable for precise navigation. Adapting LEO communication satellites for GNSS would require introducing dedicated navigation signal structures, managing interference between navigation and communication channels, and potentially reallocating onboard power resources — all of which complicate satellite design and operation.

## 4. Time Synchronization and Clock Stability

**(a) Frequent Clock Corrections**

Accurate GNSS positioning depends on precise time synchronization at the nanosecond level. Traditional GNSS satellites utilize highly stable atomic clocks requiring relatively infrequent updates. In contrast, LEO satellites experience higher orbital perturbations due to atmospheric drag and solar radiation pressure, demanding more frequent clock corrections, possibly on an hourly basis, thereby increasing ground segment complexity.

**(b) Dependence on Inter-Satellite Links**

Maintaining synchronization across dense LEO constellations may require real-time inter-satellite links (ISLs), such as optical or radio frequency links. Although some communication constellations (e.g., Starlink) already employ ISLs, adapting them to meet strict GNSS timing requirements introduces additional engineering challenges, especially regarding link stability and precision.

## 5. Atmospheric and Environmental Impacts

**(a) Ionospheric Delay Variations**

Although the ionospheric path length for LEO signals is shorter than that for MEO signals, the rapid relative motion of LEO satellites leads to fast-changing signal paths and ionospheric conditions. Real-time, dynamic ionospheric correction models are thus required to maintain high positioning accuracy, increasing system complexity.

**(b) Multipath and Signal Blockage in Urban Environments**

Due to lower elevation angles, LEO satellite signals are more susceptible to multipath effects and signal blockage caused by urban structures. This can severely degrade navigation performance in critical environments such as cities, where robust and accurate positioning is most needed.

## 6. Orbital Maintenance and Space Debris Risks

**(a) Atmospheric Drag and Orbital Decay**

LEO satellites are continuously affected by atmospheric drag, leading to gradual orbital decay if not actively counteracted by propulsion systems. Frequent orbit maintenance is necessary to preserve the integrity of the navigation constellation and avoid position prediction errors.

**(b) Space Debris and Collision Risks**

The LEO environment is increasingly congested with operational satellites and debris fragments. Collisions pose serious risks to satellite survivability and navigation service continuity. Active debris monitoring, avoidance maneuvers, and post-collision recovery strategies are required, further increasing operational costs and complexity.

## 7. Regulatory and Standardization Challenges

**(a) Frequency Allocation Conflicts**

LEO communication satellites typically operate in Ku and Ka bands, which are not traditionally assigned for navigation. Allocating new frequencies for navigation use in LEO requires international coordination through organizations like the International Telecommunication Union (ITU), and may face resistance from existing spectrum holders.

**(b) National Security and Political Considerations**

Navigation services have strategic importance for military and civil operations. Some countries may impose restrictions on foreign-operated LEO-based navigation systems due to national security concerns. Consequently, achieving global acceptance would require diplomatic efforts and adherence to various national and international regulations.


The following table categorizes and summarizes the principal challenges of employing LEO communication satellites for navigation, highlighting their corresponding technical impacts.

## Summary Table of Challenges

| Category | Challenges | Key Implications |
|:---------|:-----------|:-----------------|
| **Orbital Dynamics** | Rapid satellite movement; short visibility periods | Frequent handovers; increased receiver complexity |
| **Coverage** | Small coverage footprint; massive constellation size required | High deployment and maintenance costs |
| **Signal Processing** | Severe Doppler shifts; non-navigation signals | Need for advanced receiver design and signal redesign |
| **Time Synchronization** | Frequent clock corrections; reliance on ISLs | Increased operational complexity |
| **Environmental Impacts** | Ionospheric variability; urban multipath | Degraded accuracy in dynamic or obstructed environments |
| **Orbital Maintenance** | Atmospheric drag; space debris collision risk | Requirement for active orbit management and debris avoidance |
| **Regulatory Issues** | Spectrum allocation conflicts; political sensitivities | Challenges in international adoption and coordination |

## 8. Conclusion

While LEO communication satellites present exciting opportunities for enhancing GNSS services — particularly in urban environments, autonomous navigation, and emergency response — significant challenges must be addressed before they can be reliably used for navigation. These challenges span orbital dynamics, signal design, synchronization complexity, environmental interference, collision risks, and regulatory constraints. Future advancements are expected to focus on hybrid constellation architectures (combining LEO, MEO, and GEO satellites), robust Doppler compensation techniques, improved onboard clock technologies, and international spectrum coordination. If these barriers are successfully overcome, LEO-based navigation could serve as a powerful complement to existing GNSS infrastructure, offering enhanced resilience, redundancy, and positioning performance.

---

# Task 5 – GNSS in Remote Sensing: A Focus on GNSS Interferometric Reflectometry (GNSS-IR)

## 1. Introduction

Global Navigation Satellite Systems (GNSS), including GPS (USA), Galileo (EU), and BeiDou (China), are primarily associated with positioning, navigation, and timing (PNT) services. However, recent advancements in signal processing have enabled GNSS signals to be repurposed for remote sensing applications. Among these, GNSS Interferometric Reflectometry (GNSS-IR) has emerged as a promising technique for cost-effective, continuous, and weather-resilient environmental monitoring. This essay focuses on the principles, applications, impacts, challenges and limitations, and future prospects of GNSS-IR within the context of remote sensing.

## 2. Principles of GNSS-IR

GNSS-IR leverages multipath interference caused by GNSS signals reflecting off Earth's surfaces and interfering with direct signals at the receiver.

- Signal Interference Mechanism:

   Receivers capture both direct and reflected signals. The resulting interference patterns in the Signal-to-Noise Ratio (SNR) data encode information about the reflecting surface.

- Key Observables:

  - The frequency of SNR oscillations relates to the vertical distance between the antenna and the reflecting surface, enabling the measurement of snow depth, water levels, and surface elevation changes.

  - The amplitude of SNR oscillations corresponds to surface reflectivity, which is sensitive to soil moisture, vegetation density, and surface roughness.

Through these observables, GNSS-IR passively retrieves valuable geophysical parameters critical to environmental remote sensing.

## 3. Applications of GNSS-IR in Remote Sensing

GNSS-IR supports diverse remote sensing applications:

- **Soil Moisture Monitoring**: Variations in the soil's dielectric constant affect reflected GNSS signals. GNSS-IR enables continuous estimation of surface soil moisture, supporting agricultural management and drought monitoring.
- **Snow Depth and Ice Surface Monitoring**: Changes in snow accumulation modify the antenna-to-surface height, detectable via SNR pattern shifts. GNSS-IR is particularly valuable in polar and alpine regions.
- **Water Level Monitoring**: Reflections from lakes, rivers, and coastal surfaces allow GNSS-IR to monitor water level variations as a complement to traditional tide gauges and altimetry satellites.
- **Vegetation and Biomass Assessment**: Growing vegetation affects surface reflectivity. GNSS-IR can indirectly monitor biomass and crop dynamics, providing supplementary information to optical indices such as NDVI.

These applications highlight GNSS-IR’s unique capacity to provide continuous, ground-based environmental observations.

## 4. Impact of GNSS-IR on Remote Sensing

GNSS-IR has had a notable impact on the field of remote sensing by:

- **Enhancing Temporal Resolution**: Continuous GNSS signal availability supports near-real-time monitoring of environmental processes, overcoming the revisit time limitations of traditional satellites.
- **Reducing Observation Costs**: By utilizing existing GNSS infrastructure, GNSS-IR minimizes the need for dedicated sensing platforms, democratizing access to environmental data collection.
- **Improving Observation Robustness**: GNSS signals are resilient to atmospheric conditions, enabling data acquisition under cloud cover, rain, and snow, unlike optical and infrared systems.
- **Facilitating Multi-Parameter Monitoring**: A single GNSS-IR station can simultaneously monitor soil moisture, snow depth, water levels, and vegetation growth, offering comprehensive environmental insights.

These impacts demonstrate how GNSS-IR broadens the capabilities and accessibility of terrestrial remote sensing.

## 5. Challenges and Limitations in Remote Sensing Applications

Despite its advantages, GNSS-IR faces specific challenges and inherent limitations within remote sensing:

- **Localized Spatial Coverage**: Observations are restricted to a small area surrounding each receiver (typically tens to hundreds of meters), limiting large-scale spatial representation.
- **Environmental Sensitivity**: Surface roughness, vegetation density, and terrain complexity can distort reflection signals, affecting retrieval accuracy.
- **Signal Contamination in Urban Areas**: Multipath reflections from buildings and infrastructure introduce noise, complicating the extraction of useful environmental signals.
- **Calibration and Validation Needs**: Accurate parameter retrieval often depends on ground-truth datasets, especially for soil moisture and snow depth estimates, increasing operational complexity.
- **Dependence on Surface Conditions**: Dynamic changes in surface wetness, vegetation, and seasonal snow impact signal quality, requiring adaptive processing strategies.

These factors must be addressed to fully realize the operational potential of GNSS-IR for remote sensing purposes.

## 6. Future Prospects for GNSS-IR in Remote Sensing

The future development of GNSS-IR for remote sensing is promising:

- **Expansion of GNSS Networks**: The deployment of additional GNSS-IR-capable stations will improve spatial coverage, enhancing environmental monitoring networks.
- **Integration with AI**: Machine learning algorithms will automate signal processing, noise reduction, and parameter inversion, boosting efficiency and accuracy.
- **Multi-Sensor Fusion**: Combining GNSS-IR with InSAR, optical remote sensing, and LiDAR will enable more comprehensive Earth system monitoring.
- **Utilization of Emerging GNSS Signals**: The increasing availability of multi-frequency and multi-constellation GNSS signals (e.g., BeiDou-3, Galileo) will enhance data redundancy and improve retrieval robustness.

Through these advancements, GNSS-IR is expected to become a critical component of next-generation environmental sensing systems.

## 7. Conclusion

GNSS Interferometric Reflectometry demonstrates the innovative adaptation of GNSS signals for remote sensing purposes. By exploiting multipath interference patterns, GNSS-IR enables low-cost, continuous, and weather-resilient monitoring of environmental parameters such as soil moisture, snow depth, water levels, and vegetation dynamics. Although limitations related to spatial resolution and surface variability persist, ongoing technological advancements are positioning GNSS-IR as a critical component of future Earth observation systems.

---

# References

[1] Enge, P.K. The Global Positioning System: Signals, measurements, and performance. Int J Wireless Inf Networks 1, 83–105 (1994). [https://doi.org/10.1007/BF02106512](https://doi.org/10.1007/BF02106512).

[2] Elliott Kaplan; Christopher Hegarty, Understanding GPS/GNSS: Principles and Applications, Third Edition , Artech, 2017. [https://ieeexplore.ieee.org/document/9100468](https://ieeexplore.ieee.org/document/9100468).


[3] Nop, N. DIFFERENTIAL GNSS (GLOBAL NAVIGATION SATELLITE SYSTEMS) SERVICES – VERSION 3. [https://www.ge0mlib.com/papers/Protocols/RTCM_SC-104_v3.1.pdf](https://www.ge0mlib.com/papers/Protocols/RTCM_SC-104_v3.1.pdf).

[4] Zumberge, J. F., M. B. Heflin, D. C. Jefferson, M. M. Watkins, and F. H. Webb (1997), Precise point positioning for the efficient and robust analysis of GPS data from large networks, J. Geophys. Res., 102(B3), 5005–5017, [doi:10.1029/96JB03860](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/96JB03860).

[5] Li, X., Huang, J., Li, X. et al. Review of PPP–RTK: achievements, challenges, and opportunities. Satell Navig 3, 28 (2022). [https://doi.org/10.1186/s43020-022-00089-9](https://doi.org/10.1186/s43020-022-00089-9)

[6] J. Morton, F. van Diggelen, S. Lo, and K. Pesyna, *Position, Navigation, and Timing Technologies in the 21st Century*, Wiley-IEEE Press, 2020. [https://onlinelibrary.wiley.com/doi/book/10.1002/9781119458449](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119458449).

[7] Walter, T., & Enge, P. (1995, September). Weighted RAIM for precision approach. In Proceedings of Ion GPS (Vol. 8, No. 1, pp. 1995-2004). Institute of Navigation. [https://web.stanford.edu/group/scpnt/gpslab/pubs/papers/Walter_IONGPS_1995_wraim.pdf](https://web.stanford.edu/group/scpnt/gpslab/pubs/papers/Walter_IONGPS_1995_wraim.pdf).

[8] Lei Wang, Deren Li, Ruizhi Chen, Wenju Fu, Xin Shen, Jiang Hao 2,. Low Earth Orbiter (LEO) Navigation Augmentation: Opportunities and Challenges. Strategic Study of CAE, 2020, 22(2): 144‒152. [https://doi.org/10.15302/J-SSCAE-2020.02.018](https://doi.org/10.15302/J-SSCAE-2020.02.018).

[9] M. Hui et al., "A Review of LEO Satellite Communication Payloads for Integrated Communication, Navigation, and Remote Sensing: Opportunities, Challenges, Future Directions," in IEEE Internet of Things Journal, [doi: 10.1109/JIOT.2025.3553942](https://ieeexplore.ieee.org/abstract/document/10945753).

[10] Yuehao Teng, Xiaolin Jia, Ge Peng, LEO navigation augmentation constellation design and precise point positioning performance analysis based on BDS-3, Advances in Space Research, Volume 72, Issue 6, 2023, Pages 1944-1960, ISSN 0273-1177, [https://doi.org/10.1016/j.asr.2023.05.018](https://www.sciencedirect.com/science/article/pii/S0273117723003745).


[11] Larson, K. M., E. E. Small, E. D. Gutmann, A. L. Bilich, J. J. Braun, and V. U. Zavorotny (2008), Use of GPS receivers as a soil moisture network for water cycle studies, Geophys. Res. Lett., 35, L24405, [doi:10.1029/2008GL036013](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008GL036013).

[12] C. C. Chew, E. E. Small, K. M. Larson and V. U. Zavorotny, "Vegetation Sensing Using GPS-Interferometric Reflectometry: Theoretical Effects of Canopy Parameters on Signal-to-Noise Ratio Data," in IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 5, pp. 2755-2764, May 2015, [doi: 10.1109/TGRS.2014.2364513](https://ieeexplore.ieee.org/document/6954462/).

[13] Nicolas Roussel, Guillaume Ramillien, Frédéric Frappart, José Darrozes, Adrien Gay, Richard Biancale, Nicolas Striebig, Vincent Hanquiez, Xavier Bertin, Damien Allain, Sea level monitoring and sea state estimate using a single geodetic receiver,
Remote Sensing of Environment, Volume 171, 2015, Pages 261-277, ISSN 0034-4257, [https://doi.org/10.1016/j.rse.2015.10.011](https://www.sciencedirect.com/science/article/pii/S0034425715301620).

[14] [GPT-4o](https://chatgpt.com/?model=gpt-4o)

[15] [DeepSeek](https://www.deepseek.com/)

