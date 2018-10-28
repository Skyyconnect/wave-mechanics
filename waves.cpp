#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <cmath>




class Vector3D
{
    public:
        float x,y,z, magnitude, frequency, iRot, jRot, kRot;
        int direction;
        std::vector <float> values;


    Vector3D(float x_, float y_, float z_, int dir_)
    {
        x = x_;
        y = y_;
        z = z_;
        direction = dir_;
    }


    float dotProductAngle(Vector3D vector, float theta)
    {
        return magnitude*vector.magnitude* cos(theta);
    }

    float dotProduct(Vector3D vector)
    {
        return (x*vector.x)+ (y*vector.y)+ (z*vector.z);
    }

    void twoDRotate(float theta, int i, int j)
    {
        float iRot = values[i]*cos(theta) - values[j]*sin(theta);
        float jRot = values[i]*sin(theta) + values[j]*cos(theta);
    }

   void threeD_Rotate_X(float theta, int i, int j, int k)
    {

        float jRot = values[i]*cos(theta) - values[j]*sin(theta);
        float kRot = values[i]*sin(theta) + values[j]*cos(theta);
    

    }

    void threeD_Rotate_Y(float theta, int i, int j, int k)
    {
        float iRot = values[i]*cos(theta) - values[j]*sin(theta);
        float kRot = values[i]*sin(theta) + values[j]*cos(theta);
       

    }


    void threeD_Rotate_Z(float theta, int i, int j, int k)
    {
        float iRot = values[i]*cos(theta) - values[j]*sin(theta);
        float jRot = values[i]*sin(theta) + values[j]*cos(theta);
       
    }

    bool isOrthogonal(Vector3D vector){
        return dotProduct(vector) == 0;
    }

    

};



class WaveResolver:Vector3D{
    protected:
        std::vector<float> vect;
        int currentIndex;
        float currentX, currentY, currentZ;
        std::string currentId;
        
    public:


        WaveResolver(std::vector<float> vect_, std::string currentId_, int currentIndex_):Vector3D(0.0, 0.0, 0.0, 1){
                vect = vect_;
                currentId = currentId_;
                currentIndex = currentIndex_;
                currentX = 0.0;
                currentY = 0.0; 
                currentZ = 0.0;
            }

        bool isEqual(Vector3D vectorA, Vector3D vectorB){
            return (vectorA.values.size() != vectorB.values.size());
        }


        int reverse_bits(int x,  int bits){
            int y = 0;
            for(int i = 0; i < bits; i++){
                y = (y<<1) | (x &1);
                x >>= 1;
            }
            return y;
        }

        std::vector<float> getTable(int n, bool odd){ // sin or cosine, sin is odd.
            std::vector<float> values; 
            for (int i = 0; i < n; i ++){
                int j = (i*i)%(n*2);
                if(odd){
                    values.push_back(sin(PI*j/n));
                }else{
                    values.push_back(cos(PI*j/n));
                }
            }
            return values;

        }

        int powerOfTwo(Vector3D &vector, Vector3D &vi, int  n ){
            int m = 1;
            while (m < n*2+1){
                m *= 2;
                vector.values.push_back(m);
                vi.values.push_back(m);
            }
            return m;
        }

        Vector3D process(Vector3D &vector, Vector3D &vi){ 
            int m = powerOfTwo(vector, vi, vector.values.size());
            std::vector<float> v = vector.values;
            std::vector<float> Vi_v = vi.values;
            float vTemp[m];
            float viTemp[m];
            Vector3D Vtemp(0.0,0.0,0.0,1);
            Vector3D Vitemp(0.0,0.0,0.0,1);
           

            for(int i = 0; i < v.size(); i++){
                vTemp[i] = v[i]*getTable(v.size(), false)[i]+ Vi_v[i]*getTable(v.size(), true)[i];
                viTemp[i] = -v[i]*getTable(v.size(), true)[i]+ Vi_v[i]*getTable(v.size(), false)[i];
                Vtemp.values.push_back(vTemp[i]);
                Vitemp.values.push_back(viTemp[i]);
            }
            
            for(int j = v.size(); j < m;  j++ ){
                vTemp[j] = viTemp[j] = 0.0; 
            }
             
            float newVector[m];
            float newVectori[m];
            Vector3D newVector3D(0.0,0.0,0.0,1);
            Vector3D newVector3Di(0.0,0.0,0.0,1);
            newVector[0] = getTable(v.size(), false)[0];
            newVectori[0] = getTable(v.size(), true)[0];
           
            for (int k = 1; k < v.size(); k++){
                newVector[k] = newVector[m-k] = getTable(v.size(), true)[k];
                newVectori[k] = newVectori[m-k] = getTable(v.size(), false)[k];
                newVector3D.values.push_back(newVector[k]);
                newVector3Di.values.push_back(newVectori[k]);

            }
            
            for (int l = v.size(); l <= m-v.size(); l++){
                newVector[l] = 0.0; 
                newVectori[l] = 0.0;
                newVector3D.values[l] = 0.0;
                newVector3Di.values[l] = 0.0;
            }
              
            float convoR[m];
            float convoi[m];
            std::vector<float>  convoRV;
            std::vector<float>  convoiV;
            for (int n =0; n < m; n++){
                convoRV.push_back(convoR[n]);        
                convoiV.push_back(convoi[n]);

            }
            
        
            convolutionComplex(Vtemp, Vitemp, newVector3D,newVector3Di,convoRV, convoiV);
            for(int m = 0; m < v.size(); m++){
                v[m] = convoR[m]*getTable(v.size(),false)[m] + convoi[m]*getTable(v.size(), true)[m];
                Vi_v[m] = -convoR[m] * getTable(v.size(), true)[m] + convoi[m]*getTable(v.size(), false)[m];
            }
            return vi;

        }

        std::vector<Vector3D> bitReverseAddressing(Vector3D &vector, Vector3D &vi, int levels){
            std::vector<Vector3D> tupleVector;
            for(int i = 0; i < vector.values.size(); i ++){
                int j = reverse_bits(i, levels);
                if(j > i){
                    float temp = vector.values[i];
                    vector.values[i] = vector.values[j];
                    vector.values[j] = temp;
                    temp = vi.values[i];
                    vi.values[i] = vi.values[j];
                    vi.values[j] = temp;
                }
            }
            tupleVector.push_back(vector);
            tupleVector.push_back(vi);
            return tupleVector;
        }


        
    void convolutionReal(Vector3D &vectorX, Vector3D &vectorY, std::vector<float> &out){ 
        Vector3D Vzeros(0.0,0.0,0.0,1);
        if(isEqual(vectorX, vectorY) && isEqual(vectorX,Vzeros)){
            float zeros[vectorX.values.size()];
            for(int i = 0; i < vectorX.values.size(); i ++){
                zeros[i] = 0.0;
                Vzeros.values.push_back(0.0);
                
                convolutionComplex( vectorX, Vzeros, vectorY, Vzeros, out, out);

            }
        }
        }



        void convolutionComplex(Vector3D &vectorX, Vector3D &vectorXi, Vector3D &vectorY, Vector3D &vectorYi, std::vector<float> &outR, std::vector<float> &outI){
                transform(vectorX, vectorXi);
                transform(vectorY, vectorYi);
                for(int i = 0; i < vectorX.values.size(); i ++){
                    float temp = vectorX.values[i] * vectorY.values[i]-vectorXi.values[i]*vectorYi.values[i];
                    vectorXi.values[i] = vectorXi.values[i]* vectorY.values[i] + vectorX.values[i]*vectorYi.values[i];
                    vectorX.values[i] = temp;
                }        
                inverseTransform(vectorX,vectorXi);
                for(int j = 0; j < vectorX.values.size(); j++){
                    outR[j] = vectorX.values[j]/vectorX.values.size();
                    outI[j] = vectorXi.values[j]/ vectorX.values.size();
                }
        }

        std::vector<Vector3D> radix(Vector3D &vector, Vector3D &vi){
            std::vector<Vector3D> tupleVector;
            int n = vector.values.size();
            if(isEqual(vector,vi)){
                if(n != 1){
                    int levels = -1;
                    for(int i = 0; i < 32; i++){
                        if (1 << i == n){
                            levels = i;
                        }
                    }
                    if(levels == -1)
                        return tupleVector;
                    std::vector<float> cosine = getTable(n/2, false);
                    std::vector<float> sine = getTable(n/2, true);
                    for (int j = 0; j < n/2; j++){
                        cosine[j] = cos(2*PI*j/(n/2));
                        sine[j] = sin(2*PI*j/(n/2));
                    } 

                    bitReverseAddressing(vector, vi, levels);
                    cooleyRadix(vector,vi);
                    tupleVector.push_back(vector);
                    tupleVector.push_back(vi);
                


                }   
            }

            return tupleVector;

        }


        std::vector<Vector3D> cooleyRadix(Vector3D &vector, Vector3D &vi){
            std::vector<Vector3D> tupleVector;
            for(int size = 2; size <= vector.values.size(); size *=2){
                int halfsize = size/2;
                float tableStep = vector.values.size()/ size;
                for(int i = 0; i < vector.values.size(); i++){
                    int k = 0;
                    for (int j = i; j < i+halfsize; j++){
                        k +=tableStep;
                        float tpre = vector.values[j+halfsize]*getTable(vector.values.size(), false)[k] + vi.values[j + halfsize]* getTable(vector.values.size(), true)[k];
                        float tpim = - vector.values[j+halfsize]*getTable(vector.values.size(), true)[k] + vi.values[j + halfsize]* getTable(vector.values.size(), false)[k];
                        vector.values[j+halfsize] = vector.values[j] - tpre;
                        vi.values[j+halfsize] = vi.values[j] - tpim;
                        vector.values[j] += tpre;
                        vi.values[j] += tpim;
                    }
                }
            }
            tupleVector.push_back(vector);
            tupleVector.push_back(vi);
            return tupleVector;
        }


        std::vector<Vector3D> transform(Vector3D &vector, Vector3D &vi){
            std::vector<Vector3D> tupleVector;
            if(isEqual(vector, vi)){
                if(vector.values.size() != 0){
                    if(vector.values.size() && vector.values.size()-1 == 0){
                        vector = radix(vector, vi)[0];
                        vi = radix(vector,vi)[1];

                    }else{
                        vector = process(vector, vi);
                        vi = process(vector, vi);
                    }
                }
            }
            tupleVector.push_back(vector);
            tupleVector.push_back(vi);
            return tupleVector;
        }

        void inverseTransform(Vector3D &vector, Vector3D &vi){
            transform(vi, vector);
        }



};