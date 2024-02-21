/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "Gaudi/Property.h"
#include "GaudiAlg/Producer.h"

// Define BaseClass_t
#include "k4FWCore/BaseClass.h"

#include "edm4hep/MCParticleCollection.h"

#include <string>

struct CreateExampleEventData final : Gaudi::Functional::Producer<edm4hep::MCParticleCollection(), BaseClass_t> {
  CreateExampleEventData(const std::string& name, ISvcLocator* svcLoc)
      : Producer(name, svcLoc, KeyValue("OutputLocation", "ExampleParticles")) {}

  edm4hep::MCParticleCollection operator()() const override {
    auto coll     = edm4hep::MCParticleCollection();
    auto particle = coll.create();
    particle.setMass(m_exampleEnergy);
    info() << "Producing a particle with mass = " << m_exampleEnergy << " GeV" << endmsg;
    return coll;
  }

  Gaudi::Property<double> m_exampleEnergy{this, "ExampleEnergy", 3, "Example int to be produced"};
};

DECLARE_COMPONENT(CreateExampleEventData)
