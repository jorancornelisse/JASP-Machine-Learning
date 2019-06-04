//
// Copyright (C) 2013-2018 University of Amsterdam
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//

import QtQuick 2.8
import QtQuick.Layouts 1.3
import JASP.Controls 1.0
import JASP.Widgets 1.0

Form {

    VariablesForm {
        AvailableVariablesList {name: "variables"}
        AssignedVariablesList {
            name: "predictors"
            title: qsTr("Variables")
            singleVariable: false
            allowedColumns: ["ordinal", "scale"]
        }
    }

    GroupBox {
        title: qsTr("Tables")

        CheckBox { text: qsTr("Cluster information") ; name: "tableClusterInformation" ; enabled: true ; id: clusterInfo; checked: true
          }
    }

    GroupBox {
        title: qsTr("Plots")
        CheckBox { text: qsTr("Cluster plot")       ; name: "plot2dCluster" ; checked: false; enabled: true}
        CheckBox { text: qsTr("Dendrogram")         ; name: "dendrogram"    ; checked: false; enabled: true}
    }

    Section {
        title: qsTr("Training parameters")

        GridLayout {
          RadioButtonGroup {
              title: qsTr("Model optimization")
              name: "modelOpt"
              RadioButton { text: qsTr("Manual")                          ; name: "validationManual"; id: validationManual }
              RadioButton { text: qsTr("Gap")                             ; name: "gap"   ; id: gap }
              RadioButton { text: qsTr("Silhouette")                      ; name: "silhouette"  ; id: silhouette }
          }

          GroupBox {
              IntegerField { name: "noOfClusters"; text: qsTr("Clusters:") ; defaultValue: 3 ; min: 1; max: 999999; fieldWidth: 60; enabled: validationManual.checked }
              IntegerField { name: "maxClusters"; text: qsTr("Max. clusters:") ; defaultValue: 10 ; min: 1; max: 999999; fieldWidth: 60; enabled: validationManual.checked ? false : true }
              CheckBox { text: qsTr("Scale variables") ; name: "scaleEqualSD"; checked: true}
              CheckBox { name: "seedBox"; text: qsTr("Set seed:"); childrenOnSameRow: true
                  DoubleField  { name: "seed"; defaultValue: 1; min: -999999; max: 999999; fieldWidth: 60 }
              }
              ComboBox { name: "distance"; label: qsTr("Distance metric:");
                  model: ListModel {
                      ListElement { key: "Euclidean"            ; value: "Euclidean" }
                      ListElement { key: "Pearson correlation"                    ; value: "Pearson correlation" }
                  }
              }
              ComboBox { name: "linkage"; label: qsTr("Linkage:");
                  model: ListModel {
                      ListElement { key: "average"                    ; value: "average" }
                      ListElement { key: "single"                     ; value: "single" }
                      ListElement { key: "complete"                   ; value: "complete" }
                      ListElement { key: "centroid"                   ; value: "centroid" }
                  }
              }
          }
        }
      }

      Section {
        text: qsTr("Predictions")

            RadioButtonGroup
            {
                name: "applyModel"
                RadioButton { value: "noApp"         ; text: qsTr("Do not predict data"); checked: true        }
                RadioButton { value: "applyImpute"   ; text: qsTr("Predict missing values in target")  }
                RadioButton { value: "applyIndicator"; text: qsTr("Predict data according to apply indicator"); id: applyIndicator       }
            }

            VariablesForm {
            visible: applyIndicator.checked
                height: 150
                AvailableVariablesList { name: "predictionVariables"; allowedColumns: ["nominal"] }
                AssignedVariablesList {
                            name: "indicator"
                            title: qsTr("Apply indicator")
                            singleVariable: true
                            allowedColumns: ["nominal"]
                        }
            }
      }
}
