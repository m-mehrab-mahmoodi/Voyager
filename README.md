# Voyager 🚀

**Voyager** is an interactive 2D physics simulation and visualization desktop application built with **C++** and **Qt 6**. It allows users to simulate particle dynamics under electromagnetic fields, record collisions, and plot real-time kinematic analytics.

---

## 🌟 Key Features

### 🔬 Physics & Simulation Engine
- **Field Interactions:** Support for customized Electric and Magnetic field dynamics.
- **Particle Kinematics:** Real-time updates for position, mass, electric charge, and velocity vectors.
- **Collision Analytics:** Logging and history tracking for particle collisions.

### 🎨 Graphical User Interface (GUI) & Analytics
- **Interactive Control Panel:** Dedicated dialogs for adding/editing objects and custom fields on the fly.
- **Real-Time Data Plotting:** Integrated 2D charting (`QChartView`) to plot speed vs. time for individual objects.
- **Custom UI Styling:** Dynamic GUI theme palette loaded via external XML configurations.
- **Simulation Controls:** Step-by-step frame controls (Forward/Step Back) with configurable time step ($dt$).

---

## 🛠️ Built With

* **Language:** C++17
* **Framework:** Qt 6 (QtWidgets, QtCharts)
* **Build System:** CMake

---

## 🚀 Getting Started

### Prerequisites
Make sure you have the following installed:
* **Qt 6.x** SDK (with `QtCharts` module)
* **CMake** (v3.16 or higher)
* C++ Compiler with C++17 support (GCC, Clang, or MSVC)

### Building from Source

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/m-mehrab-mahmoodi/Voyager.git](https://github.com/m-mehrab-mahmoodi/Voyager.git)
   cd Voyager