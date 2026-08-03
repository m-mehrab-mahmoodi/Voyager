#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QTextEdit>
#include <QListWidget>
#include <QTimer>
#include <QFileDialog>
#include <QMessageBox>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QThread>
#include <QInputDialog>
#include <QDialog>
#include <QFormLayout>
#include <QDialogButtonBox>
#include <QComboBox>
#include <QCheckBox>
#include <QLineEdit>
#include <QSpinBox>
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <sstream>
#include <Qdebug>
#include <QPalette>
#include <QFile>
#include <QXmlStreamReader>
#include <QColor>
#include <QBrush>
#include <QChartView>
#include <QLineSeries>
#include <QtCharts/QChart>
#include <QtCharts/QValueAxis>

namespace fs = std::filesystem;

QPalette loadPaletteFromXml(const QString &filePath) {
    QPalette palette;
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return palette;

    QXmlStreamReader xml(&file);

    QPalette::ColorGroup currentGroup = QPalette::Active;
    QPalette::ColorRole currentRole = QPalette::Window;

    while (!xml.atEnd()) {
        xml.readNext();
        if (xml.isStartElement()) {
            if (xml.name() == "active") currentGroup = QPalette::Active;
            else if (xml.name() == "inactive") currentGroup = QPalette::Inactive;
            else if (xml.name() == "disabled") currentGroup = QPalette::Disabled;
            else if (xml.name() == "colorrole") {
                auto role = xml.attributes().value("role").toString();
                currentRole = QPalette::Window;
                if (role == "Window") currentRole = QPalette::Window;
                else if (role == "Base") currentRole = QPalette::Base;
                else if (role == "AlternateBase") currentRole = QPalette::AlternateBase;
                else if (role == "ToolTipBase") currentRole = QPalette::ToolTipBase;
                else if (role == "ToolTipText") currentRole = QPalette::ToolTipText;
                else if (role == "Button") currentRole = QPalette::Button;
                else if (role == "BrightText") currentRole = QPalette::BrightText;
                else if (role == "HighlightedText") currentRole = QPalette::HighlightedText;
                else if (role == "Light") currentRole = QPalette::Light;
            }
            else if (xml.name() == "color") {
                int r = 0, g = 0, b = 0, a = 255;
                while (!(xml.isEndElement() && xml.name() == "color")) {
                    xml.readNext();
                    if (xml.isStartElement()) {
                        if (xml.name() == "red") r = xml.readElementText().toInt();
                        else if (xml.name() == "green") g = xml.readElementText().toInt();
                        else if (xml.name() == "blue") b = xml.readElementText().toInt();
                        else if (xml.name() == "alpha") a = xml.readElementText().toInt();
                    }
                }
                QColor color(r, g, b, a);
                QBrush brush(color);
                palette.setBrush(currentGroup, currentRole, brush);
            }
        }
    }

    file.close();
    return palette;
}

class object {
public:
    std::string name;
    double R = 0;
    double force[3] = {0, 0, 0};
    double q = 0;
    double b = 0;
    double mass = 0;
    double velocity[3] = {0, 0, 0};
    double x = 0;
    double y = 0;
    double z = 0;
    std::vector<double> forcex;
    std::vector<double> forcey;
    std::vector<double> forcez;
    bool flag_grave = false;

    object(std::string name, double R, const double force[3], double mass, double q, double b) {
        this->name = name;
        this->R = R;
        for (int i = 0; i < 3; i++) {
            this->force[i] = force[i];
        }
        this->mass = mass;
        this->q = q;
        this->b = b;
    }

    void resetForces() {
        forcex.clear();
        forcey.clear();
        forcez.clear();
    }

    ~object() {}
};

enum direct {
    dir_x=1,
    dir_y,
    dir_z
};

struct s_size
{
    double x_sta=0;
    double y_sta=0;
    double z_sta=0;
    double x_end=0;
    double y_end=0;
    double z_end=0;

    bool is_in_me(double x, double y, double z) {
        return (x >= x_sta && x <= x_end && y >= y_sta && y <= y_end && z >= z_end && z <= z_end);
    }
};

class fild {
public :
    fild (std::string name, double grav, double E, double B, s_size size, bool isall, direct fild_di) : name(name), grav(grav), E(E), B(B), size(size), isall(isall), fild_di(fild_di) {}

    std::string name;
    direct fild_di;
    s_size size;
    bool isall = false;
    double grav=0,E=0,B=0;
};

void forcegrav(double forceout[3], object obj_i ,object obj_j ) {
    const double G = 6.67408e-11;
    double dx = obj_j.x - obj_i.x;
    double dy = obj_j.y - obj_i.y;
    double dz = obj_j.z - obj_i.z;

    double r2 = dx*dx + dy*dy + dz*dz;

    if (r2 < 1e-10) {
        forceout[0] = forceout[1] = forceout[2] = 0.0;
        return;
    }
    double r = sqrt(r2);

    double force_val = G * obj_i.mass * obj_j.mass / r2;

    forceout[0] = force_val * (dx / r);
    forceout[1] = force_val * (dy / r);
    forceout[2] = force_val * (dz / r);
}

void force_E(double forceout[3], object obj_i ,object obj_j) {
    const double k = 8.9875517873681764e9;
    double dx = obj_j.x - obj_i.x;
    double dy = obj_j.y - obj_i.y;
    double dz = obj_j.z - obj_i.z;

    double r2 = dx*dx + dy*dy + dz*dz;

    if (r2 < 1e-10) {
        forceout[0] = forceout[1] = forceout[2] = 0.0;
        return;
    }
    double r = sqrt(r2);

    double force_val = k * obj_i.q * obj_j.q / r2;

    forceout[0] = force_val * (dx / r);
    forceout[1] = force_val * (dy / r);
    forceout[2] = force_val * (dz / r);
}

void force_B(double forceout[3], object obj_i ,object obj_j) {
    const double mu = 4 * M_PI * 1e-7;
    double dx = obj_j.x - obj_i.x;
    double dy = obj_j.y - obj_i.y;
    double dz = obj_j.z - obj_i.z;

    double r2 = dx*dx + dy*dy + dz*dz;

    if (r2 < 1e-10) {
        forceout[0] = forceout[1] = forceout[2] = 0.0;
        return;
    }
    double r = sqrt(r2);

    double force_val = mu * obj_i.b * obj_j.b / (4 * M_PI * r2);

    forceout[0] =(-1)* force_val * (dx / r);
    forceout[1] =(-1)* force_val * (dy / r);
    forceout[2] =(-1)* force_val * (dz / r);
}

void fild_force_grav(double forceout[3], fild& fild, object& obj) {
    forceout[0] = 0;
    forceout[1] = 0;
    forceout[2] = 0;

    if (fild.isall || fild.size.is_in_me(obj.x, obj.y, obj.z)) {
        double force_val = fild.grav * obj.mass;
        switch (fild.fild_di) {
        case dir_x:
            forceout[0] = force_val;
            break;
        case dir_y:
            forceout[1] = force_val;
            break;
        case dir_z:
            forceout[2] = force_val;
            break;
        default:
            break;
        }
    } else {}
}

void fild_force_E(double forceout[3], fild& fild, object& obj) {
    forceout[0] = 0;
    forceout[1] = 0;
    forceout[2] = 0;

    if (fild.isall || fild.size.is_in_me(obj.x, obj.y, obj.z)) {
        double force_val = fild.E * obj.q;
        switch (fild.fild_di) {
        case dir_x:
            forceout[0] = force_val;
            break;
        case dir_y:
            forceout[1] = force_val;
            break;
        case dir_z:
            forceout[2] = force_val;
            break;
        default:
            break;
        }
    } else {}
}

void fild_force_lorentz(double forceout[3], fild& fild, object& obj) {
    forceout[0] = 0;
    forceout[1] = 0;
    forceout[2] = 0;

    if (fild.isall || fild.size.is_in_me(obj.x, obj.y, obj.z)) {
        switch (fild.fild_di) {
        case dir_x:
            forceout[1]= obj.q*fild.B * obj.velocity[2];
            forceout[2]= (-1)*obj.q*fild.B * obj.velocity[1];
            break;
        case dir_y:
            forceout[2]= obj.q*fild.B * obj.velocity[0];
            forceout[0]= (-1)*obj.q*fild.B * obj.velocity[2];
            break;
        case dir_z:
            forceout[0]= obj.q*fild.B * obj.velocity[1];
            forceout[1]= (-1)*obj.q*fild.B * obj.velocity[0];
            break;
        default:
            break;
        }
    }
}

void default_force( object& obj) {
    obj.forcex.push_back(obj.force[0]);
    obj.forcey.push_back(obj.force[1]);
    obj.forcez.push_back(obj.force[2]);
}

double forcecall(std::vector<double>& Force){
    double sum=0;
    for (size_t i = 0; i < Force.size(); i++)
    {
        sum+=Force[i];
    }
    return sum;
}

void crash_history1(std::vector<std::stringstream>& crash_history, std::string name_1, std::string name_2, double curent_t, double x, double y, double z, double vx_1, double vy_1, double vz_1, double vx_2, double vy_2, double vz_2) {
    crash_history.emplace_back();

    crash_history.back() << "Crash detected in " << name_1 << " and " << name_2 << " at time " << curent_t << "s: "
                         << "Position(" << x << ", " << y << ", " << z << "), "
                         << "Velocity(" << vx_1 << ", " << vy_1 << ", " << vz_1 << "), "
                         << "Velocity(" << vx_2 << ", " << vy_2 << ", " << vz_2 << ")";
}

void startmove(bool is_live, std::vector<std::stringstream>& crash_history, std::vector<object>& OBj, double dt) {
    std::vector<std::vector<double>> old_pos(OBj.size());
    for (size_t i = 0; i < OBj.size(); i++) {
        old_pos[i] = {OBj[i].x, OBj[i].y, OBj[i].z};
    }

    for (size_t i = 0; i < OBj.size(); i++) {
        double total_fx = forcecall(OBj[i].forcex);
        double total_fy = forcecall(OBj[i].forcey);
        double total_fz = forcecall(OBj[i].forcez);

        double ax = total_fx / OBj[i].mass;
        double ay = total_fy / OBj[i].mass;
        double az = total_fz / OBj[i].mass;

        OBj[i].velocity[0] += ax * dt;
        OBj[i].velocity[1] += ay * dt;
        OBj[i].velocity[2] += az * dt;

        OBj[i].x += OBj[i].velocity[0] * dt;
        OBj[i].y += OBj[i].velocity[1] * dt;
        OBj[i].z += OBj[i].velocity[2] * dt;
    }

    for (size_t i = 0; i < OBj.size(); i++) {
        for (size_t j = i + 1; j < OBj.size(); j++) {
            double dx = OBj[i].x - OBj[j].x;
            double dy = OBj[i].y - OBj[j].y;
            double dz = OBj[i].z - OBj[j].z;
            double distance = sqrt(dx * dx + dy * dy + dz * dz);

            if (distance < (OBj[i].R + OBj[j].R)) {
                if (!is_live) {crash_history1(crash_history, OBj[i].name, OBj[j].name, dt,
                                   (OBj[i].x + OBj[j].x) / 2,
                                   (OBj[i].y + OBj[j].y) / 2,
                                   (OBj[i].z + OBj[j].z) / 2,
                                   OBj[i].velocity[0], OBj[i].velocity[1], OBj[i].velocity[2],
                                   OBj[j].velocity[0], OBj[j].velocity[1], OBj[j].velocity[2]);}

                double overlap = (OBj[i].R + OBj[j].R) - distance;

                if (distance < 1e-9) {
                    distance = 1e-9;
                }
                double nx = dx / distance;
                double ny = dy / distance;
                double nz = dz / distance;

                OBj[i].x += nx * overlap * 0.5;
                OBj[i].y += ny * overlap * 0.5;
                OBj[i].z += nz * overlap * 0.5;

                OBj[j].x -= nx * overlap * 0.5;
                OBj[j].y -= ny * overlap * 0.5;
                OBj[j].z -= nz * overlap * 0.5;

                double v1n = OBj[i].velocity[0] * nx + OBj[i].velocity[1] * ny + OBj[i].velocity[2] * nz;
                double v2n = OBj[j].velocity[0] * nx + OBj[j].velocity[1] * ny + OBj[j].velocity[2] * nz;

                double restitution = 1.0;
                double impulse = (-(1 + restitution) * (v1n - v2n)) / (1.0 / OBj[i].mass + 1.0 / OBj[j].mass);

                OBj[i].velocity[0] += (impulse * nx) / OBj[i].mass;
                OBj[i].velocity[1] += (impulse * ny) / OBj[i].mass;
                OBj[i].velocity[2] += (impulse * nz) / OBj[i].mass;

                OBj[j].velocity[0] -= (impulse * nx) / OBj[j].mass;
                OBj[j].velocity[1] -= (impulse * ny) / OBj[j].mass;
                OBj[j].velocity[2] -= (impulse * nz) / OBj[j].mass;
            }
        }
    }
}

double dt_choser(std::vector<object>& OBj) {
    double max_v_sq = 0;
    for (const auto& o : OBj) {
        double v_sq = o.velocity[0] * o.velocity[0] + o.velocity[1] * o.velocity[1] + o.velocity[2] * o.velocity[2];
        max_v_sq = std::max(max_v_sq, v_sq);
    }
    double max_v = sqrt(max_v_sq);
    if (max_v < 1e-6) return 0.001;
    double dt = std::min(0.1, 5.0 / max_v);
    return std::max(dt, 0.0001);
}
void save(std::vector<object>& obj, std::vector<fild>& fild, std::string filename, const fs::path& file_path) {
    fs::path full_path = file_path / filename;
    std::ofstream file(full_path);
    if (!file.is_open()) {} else {
        for (const auto& o : obj) {
            file << std::fixed << std::setprecision(6);
            file << o.name << "/" << o.R << "/"
                 << o.force[0] << "/" << o.force[1] << "/" << o.force[2] << "/"
                 << o.mass << "/" << o.q << "/" << o.b << "/"
                 << o.x << "/" << o.y << "/" << o.z << "/"
                 << o.velocity[0] << "/" << o.velocity[1] << "/"
                 << o.velocity[2] << "\n";
        }
        file << "Fld:\n";
        for (const auto& f : fild) {
            file << f.name << "/" << f.grav << "/"
                 << f.E << "/" << f.B << "/"
                 << f.size.x_sta << "/" << f.size.y_sta << "/"
                 << f.size.z_sta << "/"
                 << f.size.x_end << "/" << f.size.y_end << "/"
                 << f.size.z_end << "/"
                 << (f.isall ? "true" : "false") << "/"
                 << static_cast<int>(f.fild_di) << "\n";
        }
        file.close();

        fs::path new_path = file_path / filename;
        try {
            rename(full_path, new_path);
        } catch (const fs::filesystem_error& e) {
            std::cerr << "error " << e.what() << std::endl;
        }
    }
}

void loadfromfile(std::vector<object>& obj, std::vector<fild>& filds, const fs::path& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cout << "Error opening file: " << file_path << std::endl;
        return;
    }

    std::string line;
    bool reading_fields = false;
    while (getline(file, line)) {
        if (line == "Fld:") {
            reading_fields = true;
            continue;
        }
        if (line.empty()) continue;

        std::vector<std::string> word;
        std::stringstream ss(line);
        std::string token;
        while (getline(ss, token, '/')) {
            word.push_back(token);
        }

        if (!reading_fields) {
            if (word.size() >= 11) {
                double force_arr[3];
                force_arr[0] = std::stod(word[2]);
                force_arr[1] = std::stod(word[3]);
                force_arr[2] = std::stod(word[4]);

                object o(word[0], std::stod(word[1]), force_arr, std::stod(word[5]), std::stod(word[6]), std::stod(word[7]));
                o.x = std::stod(word[8]);
                o.y = std::stod(word[9]);
                o.z = std::stod(word[10]);
                if (word.size() >= 14) {
                    o.velocity[0] = std::stod(word[11]);
                    o.velocity[1] = std::stod(word[12]);
                    o.velocity[2] = std::stod(word[13]);
                }
                obj.push_back(o);
            }
        } else {
            if (word.size() >= 12) {
                s_size size;
                size.x_sta = std::stod(word[4]);
                size.y_sta = std::stod(word[5]);
                size.z_sta = std::stod(word[6]);
                size.x_end = std::stod(word[7]);
                size.y_end = std::stod(word[8]);
                size.z_end = std::stod(word[9]);
                bool isall = (word[10] == "true");
                direct dir = static_cast<direct>(std::stoi(word[11]));
                fild f(word[0], std::stod(word[1]), std::stod(word[2]), std::stod(word[3]), size, isall, dir);
                filds.push_back(f);
            }
        }
    }
    file.close();
}
void extract_data(std::vector<std::stringstream>& crash_history, std::vector<object>& obj, std::vector<fild>& filds, fs::path file_path, std::string filename) {
    std::string full_file_phth=file_path.string()+"/"+filename+".txt";
    std::ofstream file(full_file_phth);
    if (!file.is_open()) {
        return;
    }

    file << std::fixed << std::setprecision(3);
    file <<"objects data: " << "\n";
    for (const auto& o : obj) {
        file << "object name : " << o.name << "\n";
        file << "mass: " << o.mass << "\n";
        file << "Electric charge: " << o.q << "\n";
        file << "Magnetic charge: " << o.b << "\n";
        file << "radius: " << o.R << "\n\n\n";
        file << "force: \n";
        file << "F(x): " << o.force[0] << "    F(y): " << o.force[1] << "    F(z): " << o.force[2] << "\n";
        file << "location \n x: " << o.x << "    y: " << o.y << "    z: " << o.z << "\n";
        file << "velocity \n V(x): " << o.velocity[0] << "    V(y): " << o.velocity[1] << "    V(z): " << o.velocity[2] << "\n";
        file <<"-------------------------------------------------------------------\n";
    }
    file<< "fields data: " << "\n";
    for (const auto& f : filds) {
        file <<"-------------------------------------------------------------------\n";
        file << "field name: " << f.name << "\n";
        if (f.grav != 0) {
            file << "fild type: Gravity\n";
            file << "Gravity value : " << f.grav << "\n";
        } else if (f.E != 0) {
            file << "fild type: Electric\n";
            file << "Electric field value : " << f.E << "\n";
        } else if (f.B != 0) {
            file << "fild type: Magnetic\n";
            file << "Magnetic field value : " << f.B << "\n";
        }

        switch (f.fild_di) {
        case dir_x:
            file << "Direction: X\n";
            break;
        case dir_y:
            file << "Direction: Y\n";
            break;
        case dir_z:
            file << "Direction: Z\n";
            break;
        default:
            break;
        }

        if (f.isall) {
            file << "Field affects all objects\n";
        } else {
            file << "Field size: \n";
            file << "Start: (" << "x: " << f.size.x_sta<< ", y: " << f.size.y_sta << ", z: " << f.size.z_sta << ")\n";
            file << "End: (" << "x: " << f.size.x_end << ", y: " << f.size.y_end << ", z: " << f.size.z_end << ")\n";
        }
    }

    file << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    file << "Simulation data:\n";
    file << "Number of objects: " << obj.size() << "\n";
    file << "Number of fields: " << filds.size() << "\n";

    file << "Crash history:\n";
    for (const auto& ss : crash_history) {
        file << ss.str() << "\n";
    }
    file.close();
}

void force_manager(std::vector<object>& obj,std::vector<fild> fild){
    for(auto & i : obj){
        i.resetForces();
    }

    for (size_t i=0; i<obj.size(); i++){
        default_force( obj[i]);

        for (size_t j=i+1;j<obj.size();j++){
            if (obj[i].flag_grave || obj[j].flag_grave){}else{
                double forceout[3];
                forcegrav (forceout,obj[i],obj[j]);
                obj[i].forcex.push_back(forceout[0]);
                obj[i].forcey.push_back(forceout[1]);
                obj[i].forcez.push_back(forceout[2]);

                obj[j].forcex.push_back((-1)*forceout[0]);
                obj[j].forcey.push_back((-1)*forceout[1]);
                obj[j].forcez.push_back((-1)*forceout[2]);
            }

            if (obj[i].q==0||obj[j].q==0){}else{
                double forceout[3];
                force_E (forceout,obj[i],obj[j]);
                obj[i].forcex.push_back(forceout[0]);
                obj[i].forcey.push_back(forceout[1]);
                obj[i].forcez.push_back(forceout[2]);

                obj[j].forcex.push_back((-1)*forceout[0]);
                obj[j].forcey.push_back((-1)*forceout[1]);
                obj[j].forcez.push_back((-1)*forceout[2]);
            }

            if (obj[i].b==0 || obj[j].b==0){}else{
                double forceout[3];
                force_B (forceout,obj[i],obj[j]);
                obj[i].forcex.push_back(forceout[0]);
                obj[i].forcey.push_back(forceout[1]);
                obj[i].forcez.push_back(forceout[2]);

                obj[j].forcex.push_back((-1)*forceout[0]);
                obj[j].forcey.push_back((-1)*forceout[1]);
                obj[j].forcez.push_back((-1)*forceout[2]);
            }
        }

        for (auto f :fild) {
            if (f.grav!=0){
                double forceout[3];
                fild_force_grav(forceout,f,obj[i]);
                obj[i].forcex.push_back(forceout[0]);
                obj[i].forcey.push_back(forceout[1]);
                obj[i].forcez.push_back(forceout[2]);
            }else if (f.B!=0){
                double forceout[3];
                fild_force_lorentz(forceout,f,obj[i]);
                obj[i].forcex.push_back(forceout[0]);
                obj[i].forcey.push_back(forceout[1]);
                obj[i].forcez.push_back(forceout[2]);
            }else if (f.E!=0){
                double forceout[3];
                fild_force_E(forceout,f,obj[i]);
                obj[i].forcex.push_back(forceout[0]);
                obj[i].forcey.push_back(forceout[1]);
                obj[i].forcez.push_back(forceout[2]);
            }
        }
    }
}

void add_object(std::vector<object>& obj, const std::string& name, double R, const double force[3], double mass, double q, double b) {
    object new_obj(name, R, force, mass, q, b);
    obj.push_back(new_obj);
}

void add_fild(std::vector<fild>& filds, const std::string& name, double grav, double E, double B, s_size size, bool isall, direct fild_di) {
    fild new_fild(name, grav, E, B, size, isall, fild_di);
    filds.push_back(new_fild);
}

class SimulationWorker : public QObject {
    Q_OBJECT
public:
    bool islive = false;
    bool running = false;
    bool singleStep = false;

    std::vector<object> obj;
    std::vector<fild> FILD;
    std::vector<std::stringstream> crash_history;

    std::vector<std::vector<object>> history_obj;
    std::vector<double> history_time;

    std::vector<std::vector<object>> undone_history_obj;
    std::vector<double> undone_history_time;

    double total_simulation_time = 10.0;
    double curent_t = 0;
    double custom_dt = 0.01;
    bool auto_calculate_dt = true;

signals:
    void updateUI();
    void simulationFinished();
    void errorOccurred(const QString& message);
    void currentTimeUpdated(double time);
    void customDtUpdated(double dt_value);
    void stepLimitReached(const QString& message);

public slots:
    void runSimulation() {
        try {
            if (!running && !singleStep) {
                return;
            }

            double current_dt_val = custom_dt;
            if (auto_calculate_dt) {
                current_dt_val = dt_choser(obj);
                emit customDtUpdated(current_dt_val);
            }

            double dt_to_use = current_dt_val;
            if ((curent_t + dt_to_use) > total_simulation_time) {
                dt_to_use = total_simulation_time - curent_t;
                if (dt_to_use < 1e-9) dt_to_use = 0;
            }

            if (dt_to_use == 0 && curent_t >= total_simulation_time - 1e-9) {
                running = false;
                emit simulationFinished();
                return;
            }

            if (running) {
                undone_history_obj.clear();
                undone_history_time.clear();
            }


            force_manager(obj, FILD);
            startmove(islive, crash_history, obj, dt_to_use);
            curent_t += dt_to_use;

            history_obj.push_back(obj);
            history_time.push_back(curent_t);

            emit updateUI();
            emit currentTimeUpdated(curent_t);

            if (curent_t >= total_simulation_time - 1e-9) {
                running = false;
                emit simulationFinished();
            }

            singleStep = false;

        } catch (const std::exception& e) {
            emit errorOccurred(QString("Simulation error: ") + e.what());
        }
    }

    void setRunning(bool r) {
        running = r;
    }

    void stepForward(int steps = 1) {
        if (running) return;

        for (int i = 0; i < steps; ++i) {
            if (!undone_history_obj.empty()) {
                obj = undone_history_obj.back();
                curent_t = undone_history_time.back();

                history_obj.push_back(obj);
                history_time.push_back(curent_t);

                undone_history_obj.pop_back();
                undone_history_time.pop_back();

                emit updateUI();
                emit currentTimeUpdated(curent_t);
            } else {
                if (curent_t < total_simulation_time - 1e-9) {
                    singleStep = true;
                    runSimulation();
                    if (curent_t >= total_simulation_time - 1e-9) {
                        emit stepLimitReached("Simulation has reached its end time.");
                        break;
                    }
                } else {
                    emit stepLimitReached("Simulation has reached its end time.");
                    break;
                }
            }
        }
    }

    void stepBack(int steps = 1) {
        if (running) return;

        for (int i = 0; i < steps; ++i) {
            if (history_obj.size() > 1) {
                undone_history_obj.push_back(history_obj.back());
                undone_history_time.push_back(history_time.back());

                history_obj.pop_back();
                history_time.pop_back();

                obj = history_obj.back();
                curent_t = history_time.back();

                emit updateUI();
                emit currentTimeUpdated(curent_t);
            } else if (history_obj.size() == 1 && curent_t != 0) {
                undone_history_obj.push_back(history_obj.back());
                undone_history_time.push_back(history_time.back());

                history_obj.clear();
                history_time.clear();

                curent_t = 0;
                history_obj.push_back(obj);
                history_time.push_back(curent_t);

                emit updateUI();
                emit currentTimeUpdated(curent_t);
                emit stepLimitReached("Cannot step back further. Already at initial state.");
                break;
            }
            else {
                emit stepLimitReached("Cannot step back further. Already at initial state.");
                break;
            }
        }
    }

    void setTotalSimulationTime(double time) {
        total_simulation_time = time;
    }

    void setCustomDt(double dt_value) {
        custom_dt = dt_value;
    }

    void setAutoCalculateDt(bool auto_calc) {
        auto_calculate_dt = auto_calc;
    }

    void resetSimulationState() {
        curent_t = 0;
        crash_history.clear();
        history_obj.clear();
        history_time.clear();
        undone_history_obj.clear();
        undone_history_time.clear();

        history_obj.push_back(obj);
        history_time.push_back(curent_t);
        emit currentTimeUpdated(0.0);
        emit updateUI();
    }
};

class AddObjectDialog : public QDialog {
    Q_OBJECT
public:
    AddObjectDialog(QWidget *parent = nullptr) : QDialog(parent) {
        QFormLayout *form = new QFormLayout(this);

        nameEdit = new QLineEdit;
        massEdit = new QDoubleSpinBox;
        massEdit->setRange(0.001, 100000);
        massEdit->setValue(1.0);
        chargeEdit = new QDoubleSpinBox;
        chargeEdit->setRange(-10000, 10000);
        chargeEdit->setValue(0.0);
        magChargeEdit = new QDoubleSpinBox;
        magChargeEdit->setRange(-10000, 10000);
        magChargeEdit->setValue(0.0);
        radiusEdit = new QDoubleSpinBox;
        radiusEdit->setRange(0.01, 100);
        radiusEdit->setValue(1.0);

        xPosEdit = new QDoubleSpinBox;
        xPosEdit->setRange(-10000, 10000);
        xPosEdit->setValue(0.0);
        yPosEdit = new QDoubleSpinBox;
        yPosEdit->setRange(-10000, 10000);
        yPosEdit->setValue(0.0);
        zPosEdit = new QDoubleSpinBox;
        zPosEdit->setRange(-10000, 10000);
        zPosEdit->setValue(0.0);

        xVelEdit = new QDoubleSpinBox;
        xVelEdit->setRange(-1000, 1000);
        xVelEdit->setValue(0.0);
        yVelEdit = new QDoubleSpinBox;
        yVelEdit->setRange(-1000, 1000);
        yVelEdit->setValue(0.0);
        zVelEdit = new QDoubleSpinBox;
        zVelEdit->setRange(-1000, 1000);
        zVelEdit->setValue(0.0);

        xForceEdit = new QDoubleSpinBox;
        xForceEdit->setRange(-1000, 1000);
        xForceEdit->setValue(0.0);
        yForceEdit = new QDoubleSpinBox;
        yForceEdit->setRange(-1000, 1000);
        yForceEdit->setValue(0.0);
        zForceEdit = new QDoubleSpinBox;
        zForceEdit->setRange(-1000, 1000);
        zForceEdit->setValue(0.0);

        disableGravityCheckBox = new QCheckBox("Disable Gravity for this object");
        disableGravityCheckBox->setChecked(false);

        form->addRow("Name:", nameEdit);
        form->addRow("Mass (kg):", massEdit);
        form->addRow("Charge (C):", chargeEdit);
        form->addRow("Magnetic Charge:", magChargeEdit);
        form->addRow("Radius:", radiusEdit);
        form->addRow("Position X:", xPosEdit);
        form->addRow("Position Y:", yPosEdit);
        form->addRow("Position Z:", zPosEdit);
        form->addRow("Velocity X:", xVelEdit);
        form->addRow("Velocity Y:", yVelEdit);
        form->addRow("Velocity Z:", zVelEdit);
        form->addRow("Initial Force X:", xForceEdit);
        form->addRow("Initial Force Y:", yForceEdit);
        form->addRow("Initial Force Z:", zForceEdit);
        form->addRow(disableGravityCheckBox);

        QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
        connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
        connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);

        form->addRow(buttons);
    }

    bool getDisableGravity() const { return disableGravityCheckBox->isChecked(); }

    void setObjectData(const object& obj) {
        nameEdit->setText(QString::fromStdString(obj.name));
        massEdit->setValue(obj.mass);
        chargeEdit->setValue(obj.q);
        magChargeEdit->setValue(obj.b);
        radiusEdit->setValue(obj.R);
        xPosEdit->setValue(obj.x);
        yPosEdit->setValue(obj.y);
        zPosEdit->setValue(obj.z);
        xVelEdit->setValue(obj.velocity[0]);
        yVelEdit->setValue(obj.velocity[1]);
        zVelEdit->setValue(obj.velocity[2]);
        xForceEdit->setValue(obj.force[0]);
        yForceEdit->setValue(obj.force[1]);
        zForceEdit->setValue(obj.force[2]);
        disableGravityCheckBox->setChecked(obj.flag_grave);
    }

    QString getName() const { return nameEdit->text(); }
    double getMass() const { return massEdit->value(); }
    double getCharge() const { return chargeEdit->value(); }
    double getMagCharge() const { return magChargeEdit->value(); }
    double getRadius() const { return radiusEdit->value(); }
    double getXPos() const { return xPosEdit->value(); }
    double getYPos() const { return yPosEdit->value(); }
    double getZPos() const { return zPosEdit->value(); }
    double getXVel() const { return xVelEdit->value(); }
    double getYVel() const { return yVelEdit->value(); }
    double getZVel() const { return zVelEdit->value(); }
    double getXForce() const { return xForceEdit->value(); }
    double getYForce() const { return yForceEdit->value(); }
    double getZForce() const { return zForceEdit->value(); }

private:
    QLineEdit *nameEdit;
    QDoubleSpinBox *massEdit, *chargeEdit, *magChargeEdit, *radiusEdit;
    QDoubleSpinBox *xPosEdit, *yPosEdit, *zPosEdit;
    QDoubleSpinBox *xVelEdit, *yVelEdit, *zVelEdit;
    QDoubleSpinBox *xForceEdit, *yForceEdit, *zForceEdit;
    QCheckBox *disableGravityCheckBox;
};

class AddFieldDialog : public QDialog {
    Q_OBJECT
public:
    AddFieldDialog(QWidget *parent = nullptr) : QDialog(parent) {
        QFormLayout *form = new QFormLayout(this);

        nameEdit = new QLineEdit;
        typeCombo = new QComboBox;
        typeCombo->addItems({"Gravity", "Electric", "Magnetic"});
        valueEdit = new QDoubleSpinBox;
        valueEdit->setRange(-1000, 1000);
        valueEdit->setValue(-9.81);
        dirCombo = new QComboBox;
        dirCombo->addItems({"X", "Y", "Z"});

        allSpaceCheck = new QCheckBox("Affects all space");
        allSpaceCheck->setChecked(true);

        xStartEdit = new QDoubleSpinBox;
        xStartEdit->setRange(-10000, 10000);
        yStartEdit = new QDoubleSpinBox;
        yStartEdit->setRange(-10000, 10000);
        zStartEdit = new QDoubleSpinBox;
        zStartEdit->setRange(-10000, 10000);
        xEndEdit = new QDoubleSpinBox;
        xEndEdit->setRange(-10000, 10000);
        yEndEdit = new QDoubleSpinBox;
        yEndEdit->setRange(-10000, 10000);
        zEndEdit = new QDoubleSpinBox;
        zEndEdit->setRange(-10000, 10000);

        connect(allSpaceCheck, &QCheckBox::toggled, this, [=](bool checked){
            xStartEdit->setEnabled(!checked);
            yStartEdit->setEnabled(!checked);
            zStartEdit->setEnabled(!checked);
            xEndEdit->setEnabled(!checked);
            yEndEdit->setEnabled(!checked);
            zEndEdit->setEnabled(!checked);
        });
        allSpaceCheck->toggled(true);

        connect(typeCombo, &QComboBox::currentTextChanged, this, [=](const QString& text){
            valueEdit->setSuffix(text == "Gravity" ? " m/s²" :
                                     text == "Electric" ? " N/C" : " T");
        });

        form->addRow("Name:", nameEdit);
        form->addRow("Type:", typeCombo);
        form->addRow("Value:", valueEdit);
        form->addRow("Direction:", dirCombo);
        form->addRow(allSpaceCheck);
        form->addRow("X Start:", xStartEdit);
        form->addRow("Y Start:", yStartEdit);
        form->addRow("Z Start:", zStartEdit);
        form->addRow("X End:", xEndEdit);
        form->addRow("Y End:", yEndEdit);
        form->addRow("Z End:", zEndEdit);

        QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
        connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
        connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);

        form->addRow(buttons);
    }

    QString getName() const { return nameEdit->text(); }
    QString getType() const { return typeCombo->currentText(); }
    double getValue() const { return valueEdit->value(); }
    direct getDirection() const {
        switch (dirCombo->currentIndex()) {
        case 0: return dir_x;
        case 1: return dir_y;
        case 2: return dir_z;
        }
        return dir_x;
    }
    bool isAllSpace() const { return allSpaceCheck->isChecked(); }
    s_size getSize() const {
        s_size size;
        size.x_sta = xStartEdit->value();
        size.y_sta = yStartEdit->value();
        size.z_sta = zStartEdit->value();
        size.x_end = xEndEdit->value();
        size.y_end = yEndEdit->value();
        size.z_end = zEndEdit->value();
        return size;
    }

private:
    QLineEdit *nameEdit;
    QComboBox *typeCombo, *dirCombo;
    QDoubleSpinBox *valueEdit;
    QCheckBox *allSpaceCheck;
    QDoubleSpinBox *xStartEdit, *yStartEdit, *zStartEdit;
    QDoubleSpinBox *xEndEdit, *yEndEdit, *zEndEdit;
};

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    QTextEdit *infoText;
    QListWidget *objectList;
    QListWidget *fieldList;
    QListWidget *crashList;

    QChartView *chartView;

    SimulationWorker *worker;
    QThread *workerThread;
    QTimer *simulationTimer;

    QDoubleSpinBox *totalSimulationTimeSpinBox;
    QDoubleSpinBox *customDtSpinBox;
    QCheckBox *autoCalculateDtCheckBox;
    QSpinBox *stepCountSpinBox;

    QLabel *currentTimeLabel;

    QPushButton *startButton;
    QPushButton *backButton;
    QPushButton *forwardButton;
    QPushButton *saveButton;
    QPushButton *loadButton;
    QPushButton *exportButton;
    QPushButton *addObjectButton;
    QPushButton *removeObjectButton;
    QPushButton *addFieldButton;
    QPushButton *removeFieldButton;
    QPushButton *viewChartButton;

    void setupUI();
    void setupSampleSimulation();
    void setupConnections();
    void addObject();
    void addField();
    void removeSelectedObject();
    void removeSelectedField();
    void updateAllLists();
    void updateInfoText();
    void saveSimulation();
    void loadSimulation();
    void exportData();
    void startStopSimulation();
    void handleCustomDtChange(double value);
    void handleAutoCalculateDtChange(bool checked);
    void updateCurrentTimeDisplay(double time);
    void updateCustomDtDisplay(double dt_value);
    void editObject(QListWidgetItem *item);
    void handleStepForwardClicked();
    void handleStepBackClicked();
    void displaySpeedChart(int objectIndex);
    void hideChart();
    void onViewChartButtonClicked();
};

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent),
    worker(new SimulationWorker()),
    workerThread(new QThread())
{
    worker->moveToThread(workerThread);
    simulationTimer = new QTimer(this);
    setupUI();
    setupSampleSimulation();
    setupConnections();
    workerThread->start();
}

MainWindow::~MainWindow() {
    if (simulationTimer->isActive()) {
        simulationTimer->stop();
    }
    worker->setRunning(false);
    if (workerThread->isRunning()) {
        workerThread->quit();
        workerThread->wait();
    }
    delete worker;
    delete workerThread;
    delete simulationTimer;
}

void MainWindow::setupUI() {
    QWidget *centralWidget = new QWidget(this);
    QHBoxLayout *mainLayout = new QHBoxLayout(centralWidget);

    QVBoxLayout *leftLayout = new QVBoxLayout();

    QGroupBox *objectsGroup = new QGroupBox("Objects");
    QVBoxLayout *objectsLayout = new QVBoxLayout(objectsGroup);

    QHBoxLayout *objectButtonsLayout = new QHBoxLayout();
    addObjectButton = new QPushButton("Add Object");
    removeObjectButton = new QPushButton("Remove Selected");
    objectButtonsLayout->addWidget(addObjectButton);
    objectButtonsLayout->addWidget(removeObjectButton);
    objectsLayout->addLayout(objectButtonsLayout);

    objectList = new QListWidget();
    objectsLayout->addWidget(objectList);
    leftLayout->addWidget(objectsGroup);

    QGroupBox *fieldsGroup = new QGroupBox("Fields");
    QVBoxLayout *fieldsLayout = new QVBoxLayout(fieldsGroup);

    QHBoxLayout *fieldButtonsLayout = new QHBoxLayout();
    addFieldButton = new QPushButton("Add Field");
    removeFieldButton = new QPushButton("Remove Selected");
    fieldButtonsLayout->addWidget(addFieldButton);
    fieldButtonsLayout->addWidget(removeFieldButton);
    fieldsLayout->addLayout(fieldButtonsLayout);

    fieldList = new QListWidget();
    fieldsLayout->addWidget(fieldList);
    leftLayout->addWidget(fieldsGroup);

    QGroupBox *crashesGroup = new QGroupBox("Collision History");
    QVBoxLayout *crashesLayout = new QVBoxLayout(crashesGroup);
    crashList = new QListWidget();
    crashesLayout->addWidget(crashList);
    leftLayout->addWidget(crashesGroup);

    QVBoxLayout *rightLayout = new QVBoxLayout();

    chartView = new QChartView();
    chartView->setRenderHint(QPainter::Antialiasing);
    rightLayout->addWidget(chartView, 2);
    hideChart();

    infoText = new QTextEdit();
    infoText->setReadOnly(true);
    rightLayout->addWidget(infoText, 3);

    viewChartButton = new QPushButton("مشاهده نمودار سرعت شیء انتخاب شده");
    rightLayout->addWidget(viewChartButton);

    QGroupBox *timeControlsGroup = new QGroupBox("Simulation Time Controls");
    QFormLayout *timeControlsLayout = new QFormLayout(timeControlsGroup);

    totalSimulationTimeSpinBox = new QDoubleSpinBox;
    totalSimulationTimeSpinBox->setRange(1.0, 100000.0);
    totalSimulationTimeSpinBox->setValue(worker->total_simulation_time);
    timeControlsLayout->addRow("Total Simulation Time (s):", totalSimulationTimeSpinBox);

    customDtSpinBox = new QDoubleSpinBox;
    customDtSpinBox->setRange(0.000001, 10.0);
    customDtSpinBox->setDecimals(6);
    customDtSpinBox->setValue(worker->custom_dt);
    timeControlsLayout->addRow("Custom Time Step (dt):", customDtSpinBox);

    autoCalculateDtCheckBox = new QCheckBox("Auto Calculate dt");
    autoCalculateDtCheckBox->setChecked(worker->auto_calculate_dt);
    timeControlsLayout->addRow(autoCalculateDtCheckBox);

    currentTimeLabel = new QLabel("Current Time: 0.00 s");
    timeControlsLayout->addRow(currentTimeLabel);

    rightLayout->addWidget(timeControlsGroup);

    QHBoxLayout *controlLayout = new QHBoxLayout();
    startButton = new QPushButton("Start Simulation");
    backButton = new QPushButton("Step Back");
    forwardButton = new QPushButton("Step Forward");
    saveButton = new QPushButton("Save");
    loadButton = new QPushButton("Load");
    exportButton = new QPushButton("Export Data");

    stepCountSpinBox = new QSpinBox();
    stepCountSpinBox->setRange(1, 1000);
    stepCountSpinBox->setValue(1);
    stepCountSpinBox->setSuffix(" steps");

    controlLayout->addWidget(startButton);
    controlLayout->addWidget(backButton);
    controlLayout->addWidget(forwardButton);
    controlLayout->addWidget(stepCountSpinBox);
    controlLayout->addWidget(saveButton);
    controlLayout->addWidget(loadButton);
    controlLayout->addWidget(exportButton);
    rightLayout->addLayout(controlLayout, 1);

    mainLayout->addLayout(leftLayout, 2);
    mainLayout->addLayout(rightLayout, 3);

    setCentralWidget(centralWidget);
    setWindowTitle("Voyager");
    resize(1200, 800);

    customDtSpinBox->setEnabled(!autoCalculateDtCheckBox->isChecked());
}

void MainWindow::setupConnections() {
    connect(startButton, &QPushButton::clicked, this, &MainWindow::startStopSimulation);
    connect(backButton, &QPushButton::clicked, this, &MainWindow::handleStepBackClicked);
    connect(forwardButton, &QPushButton::clicked, this, &MainWindow::handleStepForwardClicked);
    connect(saveButton, &QPushButton::clicked, this, &MainWindow::saveSimulation);
    connect(loadButton, &QPushButton::clicked, this, &MainWindow::loadSimulation);
    connect(exportButton, &QPushButton::clicked, this, &MainWindow::exportData);

    connect(addObjectButton, &QPushButton::clicked, this, &MainWindow::addObject);
    connect(removeObjectButton, &QPushButton::clicked, this, &MainWindow::removeSelectedObject);
    connect(addFieldButton, &QPushButton::clicked, this, &MainWindow::addField);
    connect(removeFieldButton, &QPushButton::clicked, this, &MainWindow::removeSelectedField);

    connect(totalSimulationTimeSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), worker, &SimulationWorker::setTotalSimulationTime);
    connect(customDtSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::handleCustomDtChange);
    connect(autoCalculateDtCheckBox, &QCheckBox::toggled, this, &MainWindow::handleAutoCalculateDtChange);

    connect(worker, &SimulationWorker::updateUI, this, &MainWindow::updateAllLists);
    connect(worker, &SimulationWorker::updateUI, this, &MainWindow::updateInfoText);
    connect(worker, &SimulationWorker::simulationFinished, this, [this]() {
        simulationTimer->stop();
        worker->setRunning(false);
        startButton->setText("Start Simulation");
        forwardButton->setEnabled(true);
        backButton->setEnabled(true);
        stepCountSpinBox->setEnabled(true);
        QMessageBox::information(this, "Simulation", "Simulation completed! Click on an object to see its speed graph.");
    });
    connect(worker, &SimulationWorker::errorOccurred, this, [](const QString& msg) {
        QMessageBox::critical(nullptr, "Error", msg);
    });
    connect(worker, &SimulationWorker::currentTimeUpdated, this, &MainWindow::updateCurrentTimeDisplay);
    connect(worker, &SimulationWorker::customDtUpdated, this, &MainWindow::updateCustomDtDisplay);
    connect(worker, &SimulationWorker::stepLimitReached, this, [](const QString& msg) {
        QMessageBox::information(nullptr, "Simulation Step", msg);
    });
    connect(objectList, &QListWidget::itemDoubleClicked, this, &MainWindow::editObject);
    connect(viewChartButton, &QPushButton::clicked, this, &MainWindow::onViewChartButtonClicked);
    connect(simulationTimer, &QTimer::timeout, worker, &SimulationWorker::runSimulation);
}

void MainWindow::hideChart() {
    if (chartView) {
        chartView->hide();
    }
}

void MainWindow::onViewChartButtonClicked() {
    int currentIndex = objectList->currentRow();
    if (currentIndex >= 0 && static_cast<size_t>(currentIndex) < worker->obj.size()) {
        if (!worker->running && worker->history_time.size() > 1) {
            displaySpeedChart(currentIndex);
        } else {
            QMessageBox::information(this, "نمودار در دسترس نیست",
                                     "لطفاً شبیه‌سازی را برای چند گام زمانی اجرا کنید تا داده‌های تاریخی تولید شود، یا شبیه‌سازی را در صورت اجرا بودن، متوقف کنید.");
            hideChart();
        }
    } else {
        QMessageBox::warning(this, "شیء انتخاب نشده است", "لطفاً یک شیء از لیست انتخاب کنید تا نمودار سرعت آن را مشاهده کنید.");
        hideChart();
    }
}

void MainWindow::displaySpeedChart(int objectIndex) {
    QLineSeries *series = new QLineSeries();

    double maxSpeed = 0.0;
    for (size_t i = 0; i < worker->history_time.size(); ++i) {
        const auto& state = worker->history_obj[i][objectIndex];
        double vx = state.velocity[0];
        double vy = state.velocity[1];
        double vz = state.velocity[2];
        double speed = std::sqrt(vx*vx + vy*vy + vz*vz);

        series->append(worker->history_time[i], speed);

        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    QChart *chart = new QChart();

    chart->legend()->hide();
    chart->addSeries(series);

    QString chartTitle = QString("Speed vs. Time for %1")
                             .arg(QString::fromStdString(worker->obj[objectIndex].name));
    chart->setTitle(chartTitle);

    QValueAxis *axisX = new QValueAxis;
    axisX->setTitleText("Time (s)");
    axisX->setLabelFormat("%.2f");
    axisX->setTickCount(10);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis;
    axisY->setTitleText("Speed (m/s)");
    axisY->setLabelFormat("%.2f");
    axisY->setRange(0, maxSpeed * 1.1 == 0 ? 1.0 : maxSpeed * 1.1);
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    chartView->setChart(chart);
    chartView->show();
}

void MainWindow::handleStepForwardClicked() {
    hideChart();
    int steps = stepCountSpinBox->value();
    worker->stepForward(steps);
}

void MainWindow::handleStepBackClicked() {
    hideChart();
    int steps = stepCountSpinBox->value();
    worker->stepBack(steps);
}

void MainWindow::startStopSimulation() {
    hideChart();
    if (startButton->text() == "Start Simulation") {
        worker->setRunning(true);

        if (worker->curent_t >= worker->total_simulation_time - 1e-9 || worker->history_obj.empty() || worker->curent_t == 0) {
            if(!worker->history_obj.empty()){
                worker->obj = worker->history_obj.front();
            }
            worker->resetSimulationState();
        } else {
            worker->undone_history_obj.clear();
            worker->undone_history_time.clear();
        }

        simulationTimer->start(10);
        startButton->setText("Pause Simulation");
        forwardButton->setEnabled(false);
        backButton->setEnabled(false);
        stepCountSpinBox->setEnabled(false);
    } else {
        worker->setRunning(false);
        simulationTimer->stop();
        startButton->setText("Start Simulation");
        forwardButton->setEnabled(true);
        backButton->setEnabled(true);
        stepCountSpinBox->setEnabled(true);
    }
}

void MainWindow::handleCustomDtChange(double value) {
    worker->setCustomDt(value);
}

void MainWindow::handleAutoCalculateDtChange(bool checked) {
    worker->setAutoCalculateDt(checked);
    customDtSpinBox->setEnabled(!checked);
    qDebug() << "Auto Calculate Dt checkbox toggled. Checked state:" << checked;
}

void MainWindow::updateCurrentTimeDisplay(double time) {
    currentTimeLabel->setText(QString("Current Time: %1 s").arg(time, 0, 'f', 2));
}

void MainWindow::updateCustomDtDisplay(double dt_value) {
    if (worker->auto_calculate_dt) {
        customDtSpinBox->setValue(dt_value);
    }
}

void MainWindow::setupSampleSimulation() {
    hideChart();
    worker->obj.clear();
    worker->FILD.clear();

    double force_arr[3] = {0.0, 0.0, 0.0};


    object obj1("Object 1", 1.0, force_arr, 5.0, 0, 0);
    obj1.x = 0; obj1.y = 0; obj1.z = 0;
    obj1.velocity[0] = 0.5; obj1.velocity[1] = 0.2; obj1.velocity[2] = 0;
    worker->obj.push_back(obj1);

    object obj2("Object 2", 0.8, force_arr, 3.0, 0, 0);
    obj2.x = 3; obj2.y = 0; obj2.z = 0;
    obj2.velocity[0] = -0.3; obj2.velocity[1] =5; obj2.velocity[2] = 0;
    worker->obj.push_back(obj2);

    s_size allSpace;
    s_size magneticSize;
    magneticSize.x_sta = -10; magneticSize.x_end = 10;
    magneticSize.y_sta = -10; magneticSize.y_end = 10;
    magneticSize.z_sta = -10; magneticSize.z_end = 10;

    add_fild(worker->FILD, "Gravity", -9.81, 0.0, 0.0, allSpace, true, dir_y);
    add_fild(worker->FILD, "Magnetic Field", 0.0, 0.0, 0.5, magneticSize, false, dir_z);

    worker->resetSimulationState();

    totalSimulationTimeSpinBox->setValue(worker->total_simulation_time);
    customDtSpinBox->setValue(worker->custom_dt);
    autoCalculateDtCheckBox->setChecked(worker->auto_calculate_dt);

    updateAllLists();
    updateInfoText();
}

void MainWindow::updateAllLists() {
    objectList->clear();
    for (const auto& obj : worker->obj) {
        QString info = QString("%1\nPos: (%2, %3, %4)\nVel: (%5, %6, %7)")
        .arg(QString::fromStdString(obj.name))
            .arg(obj.x, 0, 'f', 2)
            .arg(obj.y, 0, 'f', 2)
            .arg(obj.z, 0, 'f', 2)
            .arg(obj.velocity[0], 0, 'f', 2)
            .arg(obj.velocity[1], 0, 'f', 2)
            .arg(obj.velocity[2], 0, 'f', 2);
        objectList->addItem(info);
    }

    fieldList->clear();
    for (const auto& f : worker->FILD) {
        QString type;
        double val = 0;
        if (f.grav != 0) { type = "Gravity"; val = f.grav; }
        else if (f.E != 0) { type = "Electric"; val = f.E; }
        else if (f.B != 0) { type = "Magnetic"; val = f.B; }

        QString info = QString("%1\n%2: %3")
                           .arg(QString::fromStdString(f.name))
                           .arg(type)
                           .arg(val, 0, 'f', 2);
        fieldList->addItem(info);
    }

    crashList->clear();
    for (const auto& crash : worker->crash_history) {
        crashList->addItem(QString::fromStdString(crash.str()));
    }
}

void MainWindow::addObject() {
    AddObjectDialog dialog(this);
    if (dialog.exec() == QDialog::Accepted) {
        hideChart();
        std::string name = dialog.getName().toStdString();
        double mass = dialog.getMass();
        double charge = dialog.getCharge();
        double magCharge = dialog.getMagCharge();
        double radius = dialog.getRadius();
        double x = dialog.getXPos();
        double y = dialog.getYPos();
        double z = dialog.getZPos();
        double vx = dialog.getXVel();
        double vy = dialog.getYVel();
        double vz = dialog.getZVel();
        double fx = dialog.getXForce();
        double fy = dialog.getYForce();
        double fz = dialog.getZForce();
        bool disableGravity = dialog.getDisableGravity();

        double force_arr[3] = {fx, fy, fz};
        object newObj(name, radius, force_arr, mass, charge, magCharge);
        newObj.x = x;
        newObj.y = y;
        newObj.z = z;
        newObj.velocity[0] = vx;
        newObj.velocity[1] = vy;
        newObj.velocity[2] = vz;
        newObj.flag_grave = disableGravity;

        bool collisionDetected = false;
        const std::vector<object>& current_objects = worker->obj;

        for (const auto& existing_o : current_objects) {
            double dist_x = newObj.x - existing_o.x;
            double dist_y = newObj.y - existing_o.y;
            double dist_z = newObj.z - existing_o.z;

            double distance = std::sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);

            if (distance < (newObj.R + existing_o.R)) {
                collisionDetected = true;
                QMessageBox::warning(this, "Collision Error",
                                     QString("Error: Collision detected for object '%1' at (%2, %3, %4) with existing object '%5' at (%6, %7, %8).\n"
                                             "Please change location of '%1'.")
                                         .arg(QString::fromStdString(newObj.name))
                                         .arg(newObj.x, 0, 'f', 2)
                                         .arg(newObj.y, 0, 'f', 2)
                                         .arg(newObj.z, 0, 'f', 2)
                                         .arg(QString::fromStdString(existing_o.name))
                                         .arg(existing_o.x, 0, 'f', 2)
                                         .arg(existing_o.y, 0, 'f', 2)
                                         .arg(existing_o.z, 0, 'f', 2));
                break;
            }
        }

        if (!collisionDetected) {
            worker->obj.push_back(newObj);
            worker->resetSimulationState();
        }
    }
}

void MainWindow::addField() {
    AddFieldDialog dialog(this);
    if (dialog.exec() == QDialog::Accepted) {
        hideChart();
        std::string name = dialog.getName().toStdString();
        std::string type = dialog.getType().toStdString();
        double value = dialog.getValue();
        direct dir = dialog.getDirection();
        bool allSpace = dialog.isAllSpace();
        s_size size = dialog.getSize();

        double grav = 0, E = 0, B = 0;
        if (type == "Gravity") grav = value;
        else if (type == "Electric") E = value;
        else if (type == "Magnetic") B = value;

        add_fild(worker->FILD, name, grav, E, B, size, allSpace, dir);
        worker->resetSimulationState();
    }
}

void MainWindow::removeSelectedObject() {
    int currentRow = objectList->currentRow();
    if (currentRow >= 0 && static_cast<size_t>(currentRow) < worker->obj.size()) {
        hideChart();
        worker->obj.erase(worker->obj.begin() + currentRow);
        worker->resetSimulationState();
    }
}

void MainWindow::removeSelectedField() {
    int currentRow = fieldList->currentRow();
    if (currentRow >= 0 && static_cast<size_t>(currentRow) < worker->FILD.size()) {
        hideChart();
        worker->FILD.erase(worker->FILD.begin() + currentRow);
        worker->resetSimulationState();
    }
}

void MainWindow::saveSimulation() {
    QString filePath = QFileDialog::getSaveFileName(
        this, "Save Simulation", "", "Physics Files (*.phys)");

    if (!filePath.isEmpty()) {
        if (!filePath.endsWith(".phys")) {
            filePath += ".phys";
        }

        fs::path path = filePath.toStdString();
        fs::path parentDir = path.parent_path();
        std::string filename = path.filename().string();

        save(worker->obj, worker->FILD, filename, parentDir);
        QMessageBox::information(this, "Success", "Simulation saved successfully!");
    }
}

void MainWindow::loadSimulation() {
    QString filePath = QFileDialog::getOpenFileName(
        this, "Load Simulation", "", "Physics Files (*.phys)");

    if (!filePath.isEmpty()) {
        hideChart();
        worker->obj.clear();
        worker->FILD.clear();

        fs::path path = filePath.toStdString();
        loadfromfile(worker->obj, worker->FILD, path);

        worker->resetSimulationState();

        QMessageBox::information(this, "Success", "Simulation loaded successfully!");
    }
}

void MainWindow::exportData() {
    QString dirPath = QFileDialog::getExistingDirectory(
        this, "Select Export Directory");

    if (!dirPath.isEmpty()) {
        fs::path path = dirPath.toStdString();
        extract_data(worker->crash_history, worker->obj, worker->FILD, path, "simulation_data");
        QMessageBox::information(this, "Success", "Data exported successfully!");
    }
}

void MainWindow::updateInfoText() {
    QString info = QString("Current Time: %1 s\n").arg(worker->curent_t, 0, 'f', 2);
    info += QString("Total Simulation Time: %1 s\n").arg(worker->total_simulation_time, 0, 'f', 2);
    info += QString("Time Step (dt): %1\n").arg(worker->auto_calculate_dt ? "Auto" : QString::number(worker->custom_dt, 'f', 6));
    if (worker->auto_calculate_dt) {
        info += QString("Calculated dt: %1\n").arg(customDtSpinBox->value(), 0, 'f', 6);
    }
    info += "\n";

    if (worker->obj.empty()) {
        info += "No objects in simulation.\n";
    } else {
        info += "=== Objects ===\n";
        for (const auto& o : worker->obj) {
            info += QString("Name: %1\nMass: %2 kg\nCharge: %3 C\nMagnetic Charge: %4\nRadius: %5\nPosition: (%6, %7, %8)\nVelocity: (%9, %10, %11)\nInitial Force: (%12, %13, %14)\n\n")
            .arg(QString::fromStdString(o.name))
                .arg(o.mass, 0, 'f', 2)
                .arg(o.q, 0, 'f', 2)
                .arg(o.b, 0, 'f', 2)
                .arg(o.R, 0, 'f', 2)
                .arg(o.x, 0, 'f', 2)
                .arg(o.y, 0, 'f', 2)
                .arg(o.z, 0, 'f', 2)
                .arg(o.velocity[0], 0, 'f', 2)
                .arg(o.velocity[1], 0, 'f', 2)
                .arg(o.velocity[2], 0, 'f', 2)
                .arg(o.force[0], 0, 'f', 2)
                .arg(o.force[1], 0, 'f', 2)
                .arg(o.force[2], 0, 'f', 2);
        }
    }

    if (worker->FILD.empty()) {
        info += "No fields in simulation.\n";
    } else {
        info += "\n=== Fields ===\n";
        for (const auto& f : worker->FILD) {
            QString type, value;
            if (f.grav != 0) {
                type = "Gravity";
                value = QString::number(f.grav, 'f', 2) + " m/s²";
            } else if (f.E != 0) {
                type = "Electric";
                value = QString::number(f.E, 'f', 2) + " N/C";
            } else if (f.B != 0) {
                type = "Magnetic";
                value = QString::number(f.B, 'f', 2) + " T";
            }

            QString dirStr;
            switch (f.fild_di) {
            case dir_x: dirStr = "X"; break;
            case dir_y: dirStr = "Y"; break;
            case dir_z: dirStr = "Z"; break;
            }

            info += QString("Name: %1\nType: %2\nValue: %3\nDirection: %4\n")
                        .arg(QString::fromStdString(f.name))
                        .arg(type)
                        .arg(value)
                        .arg(dirStr);
            if (f.isall) {
                info += "Affects: All Space\n\n";
            } else {
                info += QString("Affects: X:[%1,%2] Y:[%3,%4] Z:[%5,%6]\n\n")
                .arg(f.size.x_sta, 0, 'f', 2).arg(f.size.x_end, 0, 'f', 2)
                    .arg(f.size.y_sta, 0, 'f', 2).arg(f.size.y_end, 0, 'f', 2)
                    .arg(f.size.z_sta, 0, 'f', 2).arg(f.size.z_end, 0, 'f', 2);
            }
        }
    }

    infoText->setText(info);
}

void MainWindow::editObject(QListWidgetItem *item) {
    int currentRow = objectList->row(item);
    if (currentRow >= 0 && static_cast<size_t>(currentRow) < worker->obj.size()) {
        hideChart();
        object& objToEdit = worker->obj[currentRow];

        AddObjectDialog dialog(this);
        dialog.setObjectData(objToEdit);

        if (dialog.exec() == QDialog::Accepted) {
            objToEdit.name = dialog.getName().toStdString();
            objToEdit.mass = dialog.getMass();
            objToEdit.q = dialog.getCharge();
            objToEdit.b = dialog.getMagCharge();
            objToEdit.R = dialog.getRadius();
            objToEdit.x = dialog.getXPos();
            objToEdit.y = dialog.getYPos();
            objToEdit.z = dialog.getZPos();
            objToEdit.velocity[0] = dialog.getXVel();
            objToEdit.velocity[1] = dialog.getYVel();
            objToEdit.velocity[2] = dialog.getZVel();
            objToEdit.force[0] = dialog.getXForce();
            objToEdit.force[1] = dialog.getYForce();
            objToEdit.force[2] = dialog.getZForce();
            objToEdit.flag_grave = dialog.getDisableGravity();

            worker->resetSimulationState();
        }
    }
}

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);



    QPalette myPalette = loadPaletteFromXml("1palet.xml");
    app.setPalette(myPalette);


    MainWindow w;
    w.setWindowIcon(QIcon("appicon.jpg"));
    w.show();
    return app.exec();
}

#include "main.moc"
